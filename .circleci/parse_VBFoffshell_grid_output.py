#!/usr/bin/env python

__doc__ = """
Print the cross section and cross section error for each channel, from the output of ./JHUGen.py run in parallel.
The main intended use is for circleci output.  Go to the VBF offshell grids job in the generator_reference or generator_syntax jobs.
Right click on "Download the full output as a file" and copy the link address.
Then run python """+__file__+""" --url https://circleci.com/...
You can also run on a local file by running without --url, and you can change the sorting options with --sort-by.
"""

import argparse
import re

from contextlib import closing
from uncertainties import ufloat
from urllib2 import urlopen

def run(fileobject, sortkey=lambda x: x):
  f = fileobject
  result = {}
  for line in f:
    if "Total unweighted xsec" in line:
      match = re.search(r"Total unweighted xsec \(used by Vegas\): *([^ ]*) *\+/- *([^ ]*) *fb", line)
      xsec = ufloat(match.group(1), match.group(2))
      next(f)
      line = next(f)
      channel, id1, id2 = (int(_) for _ in line.split()[0:3])
      result[channel] = id1, id2, xsec, xsec.s/xsec.n

  fmtheader = "{:>8} | {:>3} | {:3} | {:17} | {:10}"
  fmt       = "{:8d} | {:3d} | {:3d} | {:7.2g} | {:10.2f}"

  header = fmtheader.format("hash idx", "id1", "id2", "xsec +/- error", "error/xsec")

  rows = []
  for k, v in result.iteritems():
    rows.append((k,)+v)
  rows.sort(key=sortkey)

  print header
  print "-"*len(header)
  for row in rows:
    print fmt.format(*row)

  print

  notthere = []
  for i in xrange(1, 165):
    if not any(_[0] == i for _ in rows):
      notthere.append(i)
  if notthere:
    print "Warning: some hash elements aren't in the output:\n"+", ".join(str(_) for _ in notthere)

def runlocal(filename, *args, **kwargs):
  with open(filename) as f:
    run(f, *args, **kwargs)

def runurl(url, *args, **kwargs):
  with closing(urlopen(url)) as f:
    run(f, *args, **kwargs)

if __name__ == "__main__":
  sortkeys = {
    "xsec": lambda row: row[3].nominal_value,
    "relativeerror": lambda row: row[4],
    "hash": lambda row: row[0],
    "flavors": lambda row: row[1:2],
  }

  p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  p.add_argument("filename_or_url")
  p.add_argument("--url", action="store_true")
  p.add_argument("--sort-by", default="hash", choices=sortkeys)
  args = p.parse_args()
  if args.url:
    runurl(args.filename_or_url, sortkey=sortkeys[args.sort_by])
  else:
    runlocal(args.filename_or_url, sortkey=sortkeys[args.sort_by])
