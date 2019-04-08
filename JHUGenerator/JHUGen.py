#!/usr/bin/env python

"""
Simple python wrapper for JHUGen. In the simplest
case, calling this with arguments will just call
JHUGen with those same arguments. Run with
--submit-to batch_type to submit the job
to a queue (slurm and condor are currently
supported; contact us if you'd like us to add
another batch system).
"""

import abc
import argparse
import contextlib
import datetime
import os
import re
import sys
import subprocess
import tempfile
import textwrap

@contextlib.contextmanager
def cd(newdir):
  """http://stackoverflow.com/a/24176022/5228524"""
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(newdir))
  try:
    yield
  finally:
    os.chdir(prevdir)


def JHUGen(*jhugenargs, **kwargs):
  hideoutput = kwargs.pop("hideoutput", False)
  if kwargs: raise TypeError("Unknown kwargs: " + ", ".join(kwargs))

  if hideoutput:
    function = subprocess.check_output
    functionkwargs = {"stderr": subprocess.STDOUT}
  else:
    function = subprocess.check_call
    functionkwargs = {}

  function(["./JHUGen"] + [str(_) for _ in jhugenargs], **functionkwargs)

class JobSubmitter(object):
  __metaclass__ = abc.ABCMeta
  def __init__(self, argv, jobname, jobtime):
    self.argv = argv
    self.jobname = jobname
    self.jobtime = jobtime

  def JHUGenargs(self, array_index=None):
    result = self.args.JHUGenarg

    if array_index is not None:
      assert self.args.array_index in (None, array_index)
    else:
      array_index = self.args.array_index

    if array_index is not None:
      if 66 <= self.Process <= 69:
        if self.VBFoffsh_run is not None: raise RuntimeError("With VBFoffsh_run set to {:d}, shouldn't have an array index".format(self.VBFoffsh_run))
        result.append("VBFoffsh_run={:d}".format(array_index))
      else:
        raise RuntimeError("Process {:d} shouldn't have an array index".format(self.Process))

    return result

  def dodryrun(self, *args, **kwargs):
    JHUGenargs = self.JHUGenargs(*args, **kwargs)
    try:
      JHUGen(*JHUGenargs + ["DryRun"], hideoutput=True)
    except subprocess.CalledProcessError as e:
      print "Invalid command line: " + " ".join(JHUGenargs)
      print e.output
      raise

  def submit(self):
    if 66 <= self.Process <= 69 and self.VBFoffsh_run is None and self.args.array_index is None:
      njobs = {
        66: 164,
        67: 164,
        68: 164,
        69: 175,
      }[self.Process]
      arrayjobs=range(1, njobs+1)
      for _ in arrayjobs:
        self.dodryrun(array_index=_)
    else:
      arrayjobs = None
      self.dodryrun()

    jobrunning = self.jobrunning
    if jobrunning:
      print "job {} is already running".format(jobrunning)
      return
    return self.submittoqueue(arrayjobs=arrayjobs)

  @abc.abstractmethod
  def submittoqueue(self, arrayjobs=None): pass
  @abc.abstractproperty
  def jobindex(self):
    """
    string that indicates the job index when submitting job arrays
    os.path.expandvars will be called on this in the job
    """
  @abc.abstractproperty
  def jobindexforoutputfile(self):
    """
    string that indicates the job index when submitting job arrays
    in the output and error file names
    os.path.expandvars will be called on this in the job
    """

  def commandline(self, isarray=False):
    result = [__file__, "--on-queue", os.path.dirname(os.path.abspath(__file__))]

    for varname in ("LD_LIBRARY_PATH",):
      try:
        result += ["--set-env-var", varname, os.environ[varname]]
      except KeyError:
        pass

    if isarray:
      result += ["--array-index", self.jobindex]

    return result + self.argv

  @property
  def args(self):
    return parse_args(self.argv)

  @property
  def VBFoffsh_run(self):
    VBFoffsh_run = None
    for arg in self.args.JHUGenarg:
      if arg.startswith("VBFoffsh_run="):
        VBFoffsh_run = int(arg.split("=", 1)[1])
    return VBFoffsh_run

  @property
  def DataFile(self):
    datafile = "data/output.lhe"
    for arg in self.args.JHUGenarg:
      if arg.startswith("DataFile="):
        datafile = arg.split("=", 1)[1]

    if datafile.endswith(".lhe"):
      datafile = datafile.replace(".lhe", "")

    if self.VBFoffsh_run is not None:
      datafile = "{}_{:03d}".format(datafile, self.VBFoffsh_run)

    return datafile

  @property
  def Process(self):
    result = 0
    for arg in self.args.JHUGenarg:
      if arg.startswith("Process="):
        result = int(arg.split("=", 1)[1])
    return result

  @property
  def ReadCSmax(self):
    result = False
    for arg in self.args.JHUGenarg:
      if arg.startswith("ReadCSmax="):
        value = arg.split("=", 1)[1]
        if value.lower().startswith("t"): result = True
        if value.lower().startswith("f"): result = False
        result = int(value)
      elif arg == "ReadCSmax": result = True
      elif arg == "NoReadCSmax": result = False
      elif "ReadCSmax" in arg.split("=")[0]: raise ValueError("Didn't know how to interpret argument with ReadCSmax: "+arg)
    return result

  @property
  def isforgrid(self):
    if self.ReadCSmax: return False
    if 66 <= self.Process <= 69: return True #because not self.ReadCSmax
    if self.args.grid: return True
    return False

  @property
  def outputfile(self):
    result = self.DataFile

    if 66 <= self.Process <= 69 and self.VBFoffsh_run is None and self.args.array_index is None:
      result += "_"+self.jobindexforoutputfile

    if self.isforgrid:
      result += ".grid"

    result += ".job.out"
    return result

  @property
  def errorfile(self):
    result = self.DataFile

    if 66 <= self.Process <= 69 and self.VBFoffsh_run is None and self.args.array_index is None:
      result += "_"+self.jobindexforoutputfile

    if self.isforgrid:
      result += ".grid"

    result += ".job.err"
    return result

  @property
  def logfile(self):
    if self.isforgrid: return self.DataFile+".grid.job.log"
    return self.DataFile+".job.log"

  @abc.abstractproperty
  def jobrunning(self): pass

class JobRunner(JobSubmitter):
  def dodryrun(self, *args, **kwargs):
    "no need to do anything"

  def submittoqueue(self, arrayjobs=None):
    if self.args.on_queue:
      cdto = self.args.on_queue
    else:
      cdto = os.path.dirname(__file__)

    for name, value in self.args.set_env_var:
      os.environ[name] = os.path.expandvars(value)

    if arrayjobs is None: arrayjobs = [self.args.array_index]

    with cd(cdto):
      for array_index in arrayjobs:
        JHUGen(*self.JHUGenargs(array_index=array_index))

  @property
  def jobrunning(self): return False

  @property
  def jobindex(self):
    assert False

  @property
  def jobindexforoutputfile(self):
    assert False

class JobSubmitterSlurm(JobSubmitter):
  def submittoqueue(self, arrayjobs=None):
    sbatchcommandline = [
      "sbatch",
      "--job-name", self.jobname,
      "--time", self.formattedjobtime,
      "--output", self.outputfile,
      "--error", self.errorfile,
    ]
    if arrayjobs is not None:
      sbatchcommandline += ["--array="+",".join(arrayjobs)]
    sbatchcommandline += self.commandline(isarray=arrayjobs is not None)
    try:
      sbatchout = subprocess.check_output(
        sbatchcommandline,
        stderr=subprocess.STDOUT,
      )
    except subprocess.CalledProcessError as e:
      print e.output,
      raise
    print sbatchout,
    with open(self.logfile, "w") as f:
      f.write(sbatchout)

  @property
  def formattedjobtime(self):
    return re.sub(" days?, ", "-", str(self.jobtime))

  @property
  def jobrunning(self):
    try:
      with open(self.logfile) as f:
        jobid = f.read().split()[-1]
    except IOError:
      return False
    if jobid in subprocess.check_output(["squeue", "--job", str(jobid)]):
      return jobid
    return False

  @property
  def jobindex(self):
    return "$SLURM_ARRAY_TASK_ID"
  @property
  def jobindexforoutputfile(self):
    return "%a"

class JobSubmitterCondor(JobSubmitter):
  @staticmethod
  def quote(commandlineoption):
    return "'" + commandlineoption.replace("'", "''").replace('"', '""') + "'"
  @property
  def formattedjobtime(self):
    return "{:d}".format(int(self.jobtime.total_seconds()))
  def submittoqueue(self, arrayjobs=None):
    commandline = self.commandline(isarray=arrayjobs is not None)
    with tempfile.NamedTemporaryFile(bufsize=0) as f:
      f.write(
        textwrap.dedent(
          """
            executable = {}
            arguments = {}

            output = {self.outputfile}
            error = {self.errorfile}
            log = {self.logfile}

            +MaxRuntime = {self.formattedjobtime}
            periodic_remove = JobStatus == 5

            queue {queue}
          """.format(
            commandline[0],
            '"' + " ".join(self.quote(_) for _ in commandline[1:]) + '"',
            self=self,
            queue=(
              "1" if arrayjobs is None
              else ("JobIndexInArray in " + " ".join("{:d}".format(_) for _ in arrayjobs))
            )
          )
        )
      )
      subprocess.check_call(["condor_submit", "--batch-name", self.jobname, f.name])
  @property
  def jobrunning(self):
    jobids = set()
    try:
      with open(self.logfile) as f:
        for line in f:
          match = re.search(r"^00. \(([0-9]+[.][0-9]+)", line)
          if match:
            if line.startswith("004") or line.startswith("005") or line.startswith("009"):
              jobids.discard(match.group(1))
            elif line.startswith("000"):
              jobids.add(match.group(1))
        return ", ".join(sorted(jobids))
    except IOError:
      return False
  @property
  def jobindex(self):
    return "$(JobIndexInArray)"
  @property
  def jobindexforoutputfile(self):
    return self.jobindex


def jobsubmitter(submitto, *args, **kwargs):
  return {
    None: JobRunner,
    "slurm": JobSubmitterSlurm,
    "condor": JobSubmitterCondor,
  }[submitto](*args, **kwargs)

def submitjob(submitto, *args, **kwargs):
  submitter = jobsubmitter(submitto, *args, **kwargs)
  submitter.submit()

def main(argv):
  args = parse_args(argv)
  return submitjob(args.submit_to, argv=sys.argv[1:], jobname=args.job_name, jobtime=args.job_time)

def parse_time(time):
  #http://stackoverflow.com/a/4628148/5228524
  match = re.match(r'((?P<days>\d+?)-)?((?P<hours>\d+?):)?((?P<minutes>\d+?):)?((?P<seconds>\d+?))?', time)
  if not match: raise ValueError("Bad time string {}".format(jobtime))
  kwargs = {key: int(value) for key, value in match.groupdict().iteritems() if value is not None}
  return datetime.timedelta(**kwargs)

def parse_args(*args, **kwargs):
  p = argparse.ArgumentParser()
  p.add_argument("JHUGenarg", nargs="*")
  p.add_argument("--submit-to", choices=("condor", "slurm"))
  p.add_argument("--job-name", default="JHUGen")
  p.add_argument("--job-time", type=parse_time)
  p.add_argument("--grid", action="store_true", help="set this to true if the purpose of this job is to generate a grid to be used by later jobs")

  #The following arguments are used internally when this script
  #submits itself to a queue.  They shouldn't be set manually.
  p.add_argument("--on-queue", help=argparse.SUPPRESS)
  p.add_argument("--set-env-var", help=argparse.SUPPRESS, nargs=2, action="append", default=[])
  p.add_argument("--array-index", help=argparse.SUPPRESS, type=os.path.expandvars)

  args = p.parse_args(*args, **kwargs)

  if args.on_queue: args.submit_to = None
  if args.submit_to is not None and args.job_time is None:
    p.error("Have to supply a job time when submitting to a queue")

  return args

main(sys.argv[1:])
