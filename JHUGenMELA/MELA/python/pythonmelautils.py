from collections import OrderedDict
import os
import tempfile

import ROOT

def include(filename):
  ROOT.gROOT.ProcessLine("#include <{}>".format(filename))

def compile(filename, loadMELA=True):
  ###################################
  #this is necessary on some machines
  command = "echo '\n"
  if loadMELA: command += ".x " + os.path.join(os.path.dirname(__file__), "..", "test", "loadMELA.C+") + "\n"
  if filename: command += ".L " + filename + "+\n"
  command += "' | root -l -b "
  os.system(command)
  ###################################
  if filename: ROOT.gROOT.ProcessLine(".L {}+".format(filename))


class NamedTemporaryMacro(object):
  def __init__(self):
    self.f = tempfile.NamedTemporaryFile(suffix=".C", bufsize=0)
    self.compiled = False
  def write(self, contents):
    self.f.write(contents)
    self.compiled = False
  def compile(self):
    if not self.compiled:
      compile(self.f.name)
      self.compiled = True

ROOT.gROOT.Macro(os.path.join(os.path.dirname(__file__), "..", "test", "loadMELA.C+"))
include("Mela.h")
include("TCouplingsBase.hh")

#extract all the enums from the anonymous namespace in TCouplingsBase.hh
contents = """
  #include <TCouplingsBase.hh>
  #include <TMCFM.hh>
"""
with open(os.path.join(os.path.dirname(__file__), "..", "interface", "TCouplingsBase.hh")) as f:
  tcouplingsbasecontents = f.read()
  for enumcontents in tcouplingsbasecontents.split("enum")[1:]:
    enumcontents = enumcontents.strip()
    assert enumcontents.startswith("{")
    enumcontents = enumcontents.split("{")[1].split("}")[0]
    for enumitem in enumcontents.split(","):
      enumitem = enumitem.split("=")[0].strip()
      contents += "\n  auto py_"+enumitem+" = ::"+enumitem+";"

tcouplingsbase = NamedTemporaryMacro()
tcouplingsbase.write(contents)
tcouplingsbase.compile()

class MultiDimensionalCppArray(object):
  functionfiletemplate = """
    #include <type_traits>
    {includes}
    auto {name}_getitem({cppargs}) -> std::remove_reference<decltype({cppvariable}[item])>::type {{return {cppvariable}[item];}}
    void {name}_setitem({cppargs}, std::remove_reference<decltype({cppvariable}[item])>::type value) {{{cppvariable}[item] = value;}}
  """

  uniqueids = []
  functionfiles = {}
  getitems = {}
  setitems = {}

  def __init__(self, uniqueid, cppvariable, includes, othercppargs, *dimensions):

    self.uniqueid = uniqueid
    for i in self.uniqueids:
      if i == uniqueid:
        raise ValueError("Two MultiDimensionalCppArrays can't have the same id\n{}".format(i, self.uniqueid))
    self.uniqueids.append(uniqueid)

    self.cppvariable = cppvariable
    self.dimensions = dimensions
    self.ndim = len(self.dimensions)
    if self.ndim == 0:
      raise TypeError("Can't have a 0 dimensional array!")

    if self.ndim > 1:
      self.subarrays = []
      for i in range(dimensions[0]):
          othercppargs["int index{}".format(len(dimensions))] = i
          self.subarrays.append(
                                MultiDimensionalCppArray(
                                                         "{}_{}".format(self.uniqueid, i),
                                                         "{}[index{}]".format(cppvariable, len(dimensions)),
                                                         includes,
                                                         othercppargs,
                                                         *dimensions[1:]
                                                        )
                               )
    else:
      self.othercppargs = OrderedDict(othercppargs)
      self.functionfilecontents = self.functionfiletemplate.format(
                                                                   name="NAME",
                                                                   cppvariable=self.cppvariable,
                                                                   cppargs=",".join([key for key in self.othercppargs]+["int item"]),
                                                                   includes="\n".join("#include <{}>".format(_) for _ in includes),
                                                                  )
      self.includes = includes
      self.functionfile = self.getitem = self.setitem = None

  def writecpp(self, f=None):
    if self.ndim > 1:
      for subarray in self.subarrays:
        f = subarray.writecpp(f)
      return f

    if self.functionfilecontents not in self.functionfiles:
      if f is None:
        f = NamedTemporaryMacro()
      self.functionfiles[self.functionfilecontents] = f
      f.write(self.functionfilecontents.replace("NAME", self.uniqueid))
      return f
    else:
      return self.functionfiles[self.functionfilecontents]

  def compilecpp(self, f):
    if self.ndim > 1:
      for subarray in self.subarrays:
        subarray.compilecpp(f)
      return

    if self.functionfilecontents not in self.getitems:
      f.compile()
      self.functionfiles[self.functionfilecontents] = f
      self.getitems[self.functionfilecontents] = getattr(ROOT, "{}_getitem".format(self.uniqueid))
      self.setitems[self.functionfilecontents] = getattr(ROOT, "{}_setitem".format(self.uniqueid))

    self.functionfile = self.functionfiles[self.functionfilecontents]
    self.getitem = self.getitems[self.functionfilecontents]
    self.setitem = self.setitems[self.functionfilecontents]


  def __getitem__(self, item):
    if self.ndim > 1:
      return self.subarrays[item]
    else:
      if self.getitem is None: self.compilecpp(self.writecpp())
      if item >= self.dimensions[0]:
        raise IndexError("Index {} out of range (0-{})".format(item, self.dimensions[0]))
      return self.getitem(*(list(self.othercppargs.values())+[item]))

  def __setitem__(self, item, value):
    if self.ndim > 1:
      raise TypeError("Need to specify all indices to write to the array.")
    else:
      if self.setitem is None: self.compilecpp()
      if item >= self.dimensions[0]:
        raise IndexError("Index {} out of range (0-{})".format(item, self.dimensions[0]))
      self.setitem(*(list(self.othercppargs.values())+[item, value]))

class SelfDParameter(object):
  def __init__(self, arrayname, *indices):
    self.arrayname = arrayname
    self.indices = indices

  def __get__(self, obj, objtype):
    if obj is None: return self
    array = getattr(obj, self.arrayname)
    if not self.indices: return array
    for index in self.indices[:-1]:
      array = array[index]
    return array[self.indices[-1]]

  def __set__(self, obj, val):
    if not self.indices: return setattr(obj, self.arrayname, val)
    array = getattr(obj, self.arrayname)
    for index in self.indices[:-1]:
      array = array[index]
    array[self.indices[-1]] = val

class SelfDCoupling(object):
  def __init__(self, arrayname, *indices):
    self.arrayname = arrayname
    self.indices = indices
    self.real = SelfDParameter(arrayname, *indices+(0,))
    self.imag = SelfDParameter(arrayname, *indices+(1,))

  def __get__(self, obj, objtype):
    if obj is None: return self
    #self.real does not call __get__ on real, because it's an instance variable not a class variable
    return complex(self.real.__get__(obj, objtype), self.imag.__get__(obj, objtype))

  def __set__(self, obj, val):
    self.real.__set__(obj, val.real)
    self.imag.__set__(obj, val.imag)

