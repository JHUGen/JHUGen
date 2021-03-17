#!/usr/bin/env python

"""
This script updates the reference LHE files in reference
by downloading the LHE files created by circleci in a
generator_reference job.  This is only needed if the files
are different, i.e. if the reference_check job fails.

If the build url is
https://app.circleci.com/pipelines/github/<A>/JHUGen/<B>/workflows/<C>/jobs/<D>
then you should run
./download_reference_lhes.py <A> <D>

If the repository is private, you need to provide an access token
on the command line.  The instructions for creating one are here:
https://circleci.com/docs/2.0/managing-api-tokens/#creating-a-project-api-token
It can be done by any member of the JHUGen organization.

After running the script, the LHE files and reference_info.txt
will be updated.  Commit those changes and push.  The reference_check
job should then succeed.
"""

from __future__ import print_function

import argparse, contextlib, os, re, shutil, six
from six.moves import urllib

try:
  import yaml
except ImportError:
  raise ImportError("Please pip install pyyaml")

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("username", default="JHUGen", help="circleci username of the person or organization where the tests are run from")
p.add_argument("build_number", type=int, help="build number of the generator_reference job to get the reference LHEs from")
p.add_argument("circleci_access_token", nargs="?", help="access token for circleci if the project is private")
args = p.parse_args()

here = os.path.dirname(__file__)
with open(os.path.join(here, "config.yml")) as f:
  y = yaml.safe_load(f)

filenames = y["env"]["docker"][0]["environment"]["REFERENCE_FILENAMES"].split()

url = "https://circleci.com/api/v1.1/project/github/{}/JHUGen/{}/artifacts".format(args.username, args.build_number)
request = urllib.request.Request(url)
if args.circleci_access_token is not None:
  request.add_header("Circle-Token", args.circleci_access_token)
with contextlib.closing(urllib.request.urlopen(request)) as u:
  artifacts = re.findall('https://[^"]*', six.ensure_str(u.read()))

JHUGenversion = set()

for filename in filenames:
  relevantartifact = [artifact for artifact in artifacts if os.path.basename(artifact) == filename]
  if not relevantartifact:
    raise ValueError("Didn't find {} in the artifacts".format(filename))
  if len(relevantartifact) > 1:
    raise ValueError("Found multiple versions of {} in the artifacts".format(filename))
  relevantartifact, = relevantartifact

  request = urllib.request.Request(relevantartifact)
  if args.circleci_access_token is not None:
    request.add_header("Circle-Token", args.circleci_access_token)

  print("downloading", filename)
  with contextlib.closing(urllib.request.urlopen(request)) as copyfrom, open(os.path.join(here, "reference", filename), "wb") as copyto:
    shutil.copyfileobj(copyfrom, copyto)

  with open(os.path.join(here, "reference", filename)) as f:
    JHUGenversion |= set(re.findall(r"Output from the JHUGenerator (\S*)", f.read()))

if len(JHUGenversion) > 1:
  raise ValueError("Found inconsistent JHUGen versions in the reference lhes {}".format(JHUGenversion))
JHUGenversion, = JHUGenversion

with open(os.path.join(here, "reference", "reference_info.txt"), "w") as f:
  f.write("Downloaded from circleci {}/JHUGen build number {}\n".format(args.username, args.build_number))
  f.write("JHUGen version: "+JHUGenversion+"\n")
