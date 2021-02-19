#!/usr/bin/env python

from __future__ import print_function

import argparse, contextlib, os, re, shutil, six
from six.moves import urllib

try:
  import yaml
except ImportError:
  raise ImportError("Please pip install pyyaml")

p = argparse.ArgumentParser()
p.add_argument("username", default="JHUGen")
p.add_argument("build_number", type=int)
args = p.parse_args()

here = os.path.dirname(__file__)
with open(os.path.join(here, "config.yml")) as f:
  y = yaml.safe_load(f)

filenames = y["env"]["docker"][0]["environment"]["REFERENCE_FILENAMES"].split()

url = "https://circleci.com/api/v1.1/project/github/{}/JHUGen/{}/artifacts".format(args.username, args.build_number)
with contextlib.closing(urllib.request.urlopen(url)) as u:
  artifacts = re.findall('https://[^"]*', six.ensure_str(u.read()))

for filename in filenames:
  relevantartifact = [artifact for artifact in artifacts if os.path.basename(artifact) == filename]
  if not relevantartifact:
    raise ValueError("Didn't find {} in the artifacts".format(filename))
  if len(relevantartifact) > 1:
    raise ValueError("Found multiple versions of {} in the artifacts".format(filename))
  relevantartifact, = relevantartifact

  print("downloading", filename)
  with contextlib.closing(urllib.request.urlopen(relevantartifact)) as copyfrom, open(os.path.join(here, "reference", filename), "wb") as copyto:
    shutil.copyfileobj(copyfrom, copyto)

with open(os.path.join(here, "reference", "reference_info.txt"), "w") as f:
  f.write("Downloaded from circleci {}/JHUGen build number {}\n".format(args.username, args.build_number))
