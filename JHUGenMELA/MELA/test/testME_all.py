#!/usr/bin/env python

from collections import defaultdict
import os
import re

import ROOT

ROOT.gROOT.Macro("loadMELA.C")
ROOT.gROOT.LoadMacro("testME_v2.c+")

def defaultdictset():
  return defaultdict(set)

ratiosMCFMJHUGenVBF = defaultdict(set)
ratiosMCFMJHUGenLepZH = defaultdict(set)
ratiosMCFMJHUGenLepWH = defaultdict(set)
ratiosMCFMJHUGenHadZH = defaultdict(set)
ratiosMCFMJHUGenHadWH = defaultdict(set)

badbkgratios = set()
different = set()
badangles = set()
badselfD = set()
messedup = set()

referencefiles = os.listdir("reference")
referencefiles.sort(key=lambda _: "testME_Prop" not in _)

def unifycontent(content):
  """
  removes stuff from the reference file that would make
  it different for each run for no real reason.
  For example 0x1234567 --> 0xPOINTER
  """
  content = re.sub("0x[0-9a-f]+", "0xPOINTER", content)
#  if "Start Mela constructor" in content: content = content.split("End Mela constructor")[1].lstrip()
  return content

melaptr = ROOT.makemelaptr(13, 125, ROOT.TVar.ERROR)

for j, ref in enumerate(referencefiles, start=1):
  match = re.match("(testME_.*_Ping)_?(.*).ref", ref)
  if not match:
    raise RuntimeError("Unknown filename {}".format(ref))
  functionname = match.group(1)
  arguments = match.group(2)
  arguments = arguments.replace("2l2l", "0")
  arguments = arguments.replace("4l", "1")
  if arguments:
    arguments = [int(_.replace("TeV", "")) for _ in arguments.split("_")] + [melaptr]
  else:
    arguments = [melaptr]

  if functionname == "testME_ProdDec_MCFM_JHUGen_WBFZZ_Comparison_Ping":
    functionname = "testME_ProdDec_MCFM_JHUGen_WBFZZWW_Comparison_Ping"
    arguments.insert(3, 1)
  if functionname == "testME_ProdDec_MCFM_JHUGen_WBFWW_Comparison_Ping":
    functionname = "testME_ProdDec_MCFM_JHUGen_WBFZZWW_Comparison_Ping"
    arguments.insert(3, 2)
    arguments.insert(4, 0)

  if functionname == "testME_ProdDec_MCFM_JHUGen_WBFZZ_TU_Comparison_Ping":
    functionname = "testME_ProdDec_MCFM_JHUGen_WBFZZWW_TU_Comparison_Ping"
  if functionname == "testME_ProdDec_MCFM_JHUGen_WBFWW_TU_Comparison_Ping":
    functionname = "testME_ProdDec_MCFM_JHUGen_WBFZZWW_TU_Comparison_Ping"

  if functionname == "testME_ProdDec_MCFM_JHUGen_JJQCDZZ_Comparison_Ping":
    functionname = "testME_ProdDec_MCFM_JHUGen_JJQCDZZWW_Comparison_Ping"
    arguments.insert(1, 1)
  if functionname == "testME_ProdDec_MCFM_JHUGen_JJQCDWW_Comparison_Ping":
    functionname = "testME_ProdDec_MCFM_JHUGen_JJQCDZZWW_Comparison_Ping"
    arguments.insert(1, 2)
    arguments.insert(2, 0)

  if functionname == "testME_VH_JHUGen_13TeV_Ping":
    functionname = "testME_VH_JHUGen_Ping"
    arguments.insert(0, 13)
    arguments.insert(1, 0)

  try:
    function = getattr(ROOT, functionname)
    function(*arguments)
  except:
    print "trying to call: {}(*{})".format(functionname, arguments)
    raise

  newfile = ref.replace(".ref", ".out")

  with open(newfile) as f:
    content = unifycontent(f.read())
  with open(newfile, 'w') as f:
    f.write(content)

  with open("reference/"+ref) as f:
    refcontent = unifycontent(f.read())

  if refcontent != content:
    different.add(newfile)

  if "JHUGen/MCFM Ratio" in content:
    production = content.split("\n")[0].split()[3]
    if production == "JJVBF":
      ratios = ratiosMCFMJHUGenVBF
    elif production == "Had_ZH":
      ratios = ratiosMCFMJHUGenHadZH
    elif production in "Had_WH":
      ratios = ratiosMCFMJHUGenHadWH
    elif production == "Lep_ZH":
      ratios = ratiosMCFMJHUGenLepZH
    elif production in "Lep_WH":
      ratios = ratiosMCFMJHUGenLepWH
    else:
      assert False, production

    assert "Arrays:" in content
    sections = content.split("JHUGen/MCFM Ratio")
    for section in sections[1:]:
      sectioncontent = section.split()
      for i, number in enumerate(sectioncontent):
        #go until we find something that's not a number
        try:
          ratio = float(number)
        except ValueError:
          break
        if ratio != 0:
          ratios[ratio].add(newfile)
    assert i > 1

  if "MCFM Bkg (re-sum)/Bkg Ratio" in content:
    ratiopart = content.split("MCFM Bkg (re-sum)/Bkg Ratio")[1]
    for i, number in enumerate(ratiopart.split()):
      #go until we find something that's not a number
      try:
        ratio = float(number)
      except ValueError:
        break
      if ratio != 0 and ratio != 4:
        badbkgratios.add(newfile)
    assert i > 1

  if "angles" in content:
    for line in content.split("\n"):
      match = re.search(r"TUtil (.* angles)", line)
      if match:
        tutilline = line.split(": ")[1]
        melaline = {line for line in content.split("\n") if "MELA " + match.group(1) in line}.pop().split(": ")[1]
        if tutilline != melaline:
          badangles.add(newfile)

  if "_selfD:" in content:
    pendingselfD = {}
    for line in content.split("\n"):
      if "hadronic Z-BF" in line: break
      match = re.match("(p.*)_selfD: ([0-9.e+-]*)", line)
      if match:
        if float(match.group(2)) != pendingselfD.get(match.group(1)):
          print float(match.group(2)), pendingselfD.get(match.group(1))
        if float(match.group(2)) != pendingselfD.pop(match.group(1), None):
          badselfD.add(newfile)
      else:
        match = re.match("(p.*): ([0-9.e+-]*)", line)
        if match:
          pendingselfD[match.group(1)] = float(match.group(2))

  print j, "/", len(referencefiles)

errors = []
if different:
  errors.append("The following files are not the same as the reference files:\n  "+"\n  ".join(sorted(different)))
if len(ratiosMCFMJHUGenVBF)>1:
  errors.append("There are multiple different ratios between MCFM and JHUGen for VBF:\n" + ", ".join(str(_) for _ in ratiosMCFMJHUGenVBF))
if len(ratiosMCFMJHUGenHadZH)>1:
  errors.append("There are multiple different ratios between MCFM and JHUGen for hadronic ZH:\n" + ", ".join(str(_) for _ in ratiosMCFMJHUGenHadZH))
if len(ratiosMCFMJHUGenHadWH)>1:
  errors.append("There are multiple different ratios between MCFM and JHUGen for hadronic WH:\n" + ", ".join(str(_) for _ in ratiosMCFMJHUGenHadWH))
if len(ratiosMCFMJHUGenLepZH)>1:
  errors.append("There are multiple different ratios between MCFM and JHUGen for leptonic ZH:\n" + ", ".join(str(_) for _ in ratiosMCFMJHUGenLepZH))
if len(ratiosMCFMJHUGenLepWH)>1:
  errors.append("There are multiple different ratios between MCFM and JHUGen for leptonic WH:\n" + ", ".join(str(_) for _ in ratiosMCFMJHUGenLepWH))
if badbkgratios:
  errors.append("The following files have bad ratios between manual and automatic background:\n  "+"\n  ".join(sorted(badbkgratios)))
if badangles:
  errors.append("The following files have inconsistent angles between MELA and TUtil:\n  "+"\n  ".join(sorted(badangles)))
if badselfD:
  errors.append("The following files have inconsistent MEs between predefined and selfD hypotheses:\n  "+"\n  ".join(sorted(badselfD)))

if errors:
  raise RuntimeError("\n\n\n".join(errors))
