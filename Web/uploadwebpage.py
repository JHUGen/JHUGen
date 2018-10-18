#!/usr/bin/env python
#Inputs are down at the bottom.  No need to change anything above that unless there's a bug

import argparse
from contextlib import closing, contextmanager
from getpass import getpass
import os
import re
import shutil
from subprocess import check_call
import tempfile
import urllib2

gitdir = None

Webdir = os.getcwd()
WebGeneratordir = os.path.join(Webdir, "Generator")
if not os.path.isdir(WebGeneratordir):
     raise OSerror("You should run this from the Web directory")

def uploadwebpage(dontupload=[], dontuploadifexists=[], dryrun=False, dont_overwrite_mcfm=False, dont_overwrite_tarballs=False, edit_cadaverrc=False):
    dontupload = ListOfRegexes(dontupload)
    dontuploadifexists = ListOfRegexes(dontuploadifexists)

    if dryrun:
        JHED = password = ""
    else:
        JHED = raw_input("Enter your JHED ID: ")
        password = getpass("Enter your JHED password: ")
    if dont_overwrite_mcfm:
        dontuploadifexists.append("MCFM-precompiled")
    if dont_overwrite_tarballs:
        dontuploadifexists.append("[.]tar[.]gz$")
        dontuploadifexists.append("[.]tgz$")
    with cd("Generator"):
        create_Download(*versions)
    repmap = {
              "JHED": JHED,
              "password": password,
              "machine": "https://webdav.johnshopkins.edu/",
              "uploadfiles": getuploadfiles(".", dontupload=dontupload, dontuploadifexists=dontuploadifexists),
             }
    if dryrun:
        print "Would open cadaver, then run:"
        print
        print CadaverRC(**repmap)
    else:
        with CadaverRC(**repmap) as crc:
            if edit_cadaverrc:
                raw_input("You can now edit "+crc.filename+".  When you're done press enter.")
            with NetRC(**repmap):
                check_call(["cadaver"])


class RC(object):
    def __init__(self, **kwargs):
        self.repmap = kwargs
        if "machine" in self.repmap:
            self.repmap["machine"] = self.repmap["machine"].rstrip("/")
        self.filename
    def __enter__(self):
        if os.path.exists(self.filename): raise OSError(self.filename+" already exists!")
        with open(self.filename, "w") as f:
            f.write(self.template.format(**self.repmap))
        return self
    def __str__(self):
        return self.template.format(**self.repmap)

    def __exit__(self, exception_type, exception_value, traceback):
        os.remove(self.filename)


netrctemplate = """
default
login {JHED}
password {password}
"""

cadaverrctemplate = """
open {machine}/{JHED}/esgwebphp2/spin/public_html
ls
{uploadfiles}
exit
"""

class NetRC(RC):
    filename = os.path.expanduser("~/.netrc")
    template = netrctemplate

class CadaverRC(RC):
    filename = os.path.expanduser("~/.cadaverrc")
    template = cadaverrctemplate

def existsonwebsite(filename):
    #https://stackoverflow.com/a/16778473/5228524
    try:
        with closing(urllib2.urlopen(os.path.join(websitebase, filename))) as f:
            return True
    except urllib2.HTTPError:
        return False

class ListOfRegexes(list):
    def __contains__(self, other):
        return any(re.search(_, other) for _ in self)

def getuploadfiles(dir, dontupload, dontuploadifexists):
    files = os.listdir(dir)
    result = ""
    for f in files:
        if os.path.join(dir, f) in dontupload:
            continue
        if os.path.isfile(os.path.join(dir, f)):
            if os.path.join(dir, f) in dontuploadifexists and existsonwebsite(os.path.join(dir, f)): continue
            result += "put "+os.path.join(dir, f)+"\n"
        if os.path.isdir(os.path.join(dir, f)):
            result += "mkdir {0}\ncd {0}\n{1}cd ..\n".format(
                f,
                getuploadfiles(os.path.join(dir, f), dontupload=dontupload, dontuploadifexists=dontuploadifexists)
            )
    return result

class Version(object):
    def __init__(self, version, gitcommit=None, tarballname=None, manualname=None, manualcommit=None, melacommit=None, visible=True):
        self.version = version
        if gitcommit is None: gitcommit = version
        if tarballname is None: tarballname = "JHUGenerator.{}.tar.gz".format(self.version)
        if manualname is None: manualname = "manJHUGenerator.{}.pdf".format(self.version)
        if manualcommit is None: manualcommit = gitcommit
        self.gitcommit = gitcommit
        self.tarballname = tarballname
        self.manualname = manualname
        self.manualcommit = manualcommit
        self.melacommit = melacommit
        self.visible = visible

    def getlink(self):
        return link_template.format(tarballname=self.tarballname)

    def createtarball(self, force=False):
        global gitdir

        if not force and os.path.exists(os.path.join(WebGeneratordir, self.tarballname)): return

        if gitdir is None:
            tmpdir = tempfile.mkdtemp()
            with cd(tmpdir):
                check_call(["git", "clone", "git@github.com:JHUGen/JHUGen"])
                gitdir = os.path.join(tmpdir, "JHUGen")

        with cd(gitdir):
            check_call(["git", "clean", "-f", "-x", "-d"])
            check_call(["git", "reset", "--hard", "HEAD"])
            store = []

            #manual
            check_call(["git", "checkout", self.manualcommit])
            try:
                with cd("Manual"):
                    for i in range(3):
                        check_call(["pdflatex", "manJHUGenerator.tex"])
                    check_call(["mv", "manJHUGenerator.pdf", "../"+self.manualname])
                    store.append(self.manualname)
            except OSError: #Manual directory doesn't exist, really old version
                pass

            check_call(["git", "checkout", self.gitcommit])
            if self.melacommit is not None:
                check_call(["git", "checkout", self.melacommit, "--", "JHUGenMELA"])
            for folder in storefolders:
                if os.path.exists(folder):
                    store.append(folder)
            check_call(["tar", "-cvzf", self.tarballname] + store)
            check_call(["mv", self.tarballname, WebGeneratordir])

def create_Download(mostrecentversion, *olderversions):
    if not mostrecentversion.visible:
        raise ValueError("The most recent version has to be visible")
    mostrecentversion.createtarball(force=True)
    check_call(["cp", os.path.join(gitdir, mostrecentversion.manualname), os.path.join(Webdir, "Manual.pdf")])
    for version in [mostrecentversion] + list(olderversions) + list(reallyold):
        version.createtarball()
    Download = Download_template.format(
                                        latest = mostrecentversion.getlink(),
                                        older = "\n    ".join(v.getlink() for v in [mostrecentversion] + list(olderversions) if v.visible),
                                        reallyold = "\n    ".join(v.getlink() for v in reallyold),
                                       )
    with open(os.path.join(WebGeneratordir, "Download.html"), "w") as f:
        f.write(Download)

    #download MCFM libraries
    if not os.path.exists("MCFM-precompiled"):
        check_call(["git", "clone", "git@github.com:JHUGen/MCFM-precompiled"])
    with cd("MCFM-precompiled"):
        check_call(["git", "fetch"])
        if os.path.exists("../NNPDF30_lo_as_0130.LHgrid.tgz"):
            shutil.move("../NNPDF30_lo_as_0130.LHgrid.tgz", ".")
        check_call(["git", "checkout", MCFMprecompiledcommit])
        shutil.move("NNPDF30_lo_as_0130.LHgrid.tgz", "..")
        

@contextmanager
def cd(newdir):
    #http://stackoverflow.com/a/24176022/5228524
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

link_template = """<a href="{tarballname}">{tarballname}</a><br>"""

Download_template = """
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <head>
    <meta content="text/html; charset=ISO-8859-1"
      http-equiv="Content-Type">
    <title>Download</title>
  </head>
  <body style=" color: rgb(0, 0, 0); background-color: rgb(255, 255,
    255);" alink="#EE0000" link="#0000EE" vlink="#551A8B">
    Download:<br>
    Latest version: {latest}
    <br>
    Compiled MCFM libraries to interface with MELA that work for a good variety of software releases: <a href="MCFM-precompiled/v2/libmcfm_705.so">libmcfm_705.so</a><br>
    (See JHUGenMELA/ggZZ_MCFM/README for further information)
    <br>
    <br>
    Older generator versions:<br>
    {older}
    <br>
    Much older versions:<br>
    {reallyold}
    <br>
  </body>
</html>
""".lstrip()

########################################################################################
#inputs
#versions that should be listed on the webpage
#the top version is the most recent one
#arguments:
#   - version (mandatory), by default this is in the name of the tar ball and is also a git tag to checkout
#   the rest are optional
#   - gitcommit, what to checkout to make the tar ball.  Could be anything understood by git checkout,
#       a tag or just an SHA id (or a branch, but bad idea).  Default=version
#   - tarballname, will be used for both the link title and the name of the actual downloaded tarball.
#       Default=JHUGenerator.(version).tar.gz
#   - manualname, name of the manual pdf stored in the tarball.  Default=manJHUGenerator.(version).pdf
#   - manualcommit, name of the commit to checkout to compile the manual.  Default=gitcommit
#       could use a later version if the manual is updated after the generator is tagged
#   - melacommit, name of the commit to checkout for the JHUGenMELA folder.  Default is the same
#       as the generator
versions = (
            Version("v7.2.3"),
            Version("v7.2.2", visible=False),
            Version("v7.2.1", visible=False),
            Version("v7.2.0", manualcommit="cfb61347ffbb2391a4834f7dd803e4d533efa072", visible=False),
            Version("v7.1.4", manualcommit="e58a25479c1a16b73b657a0a6f660a27be9cd130"),
            Version("v7.1.2"),
            Version("v7.1.0", visible=False),
            Version("v7.0.11", visible=False),
            Version("v7.0.9", manualcommit="c14848e401c115e9e3f201bf7388cc750e4b09f3", melacommit="ebdb109745949008e16851674bb832c59ecc43a3", visible=False),
            Version("v7.0.2", melacommit="v6.9.8"),
            Version("v7.0.0", gitcommit="v7.0.0.beta1", manualcommit="18221e3", melacommit="v6.9.8"),
            Version("v6.9.8", manualcommit="971ad57"),
            Version("v6.9.5", manualcommit="157d32c"),
            Version("v6.8.4", gitcommit="v6.8.4.1.1"),
            Version("v6.7.8", manualcommit="3961ffa"),
            Version("v6.2.8"),
            Version("v5.6.3"),
            Version("v5.2.5", manualcommit="b9c0fa8"),
            Version("v5.2.3"),
            Version("v4.5.2"),
            Version("v4.3.2", manualcommit="4b13646"),
            Version("v4.2.1", manualcommit="41fc189"),
            Version("v4.0.1", manualcommit="757d864"),
            Version("v3.1.8"),
            Version("v3.1.1", manualcommit="0c5b754"),
            Version("v2.2.6"),
            Version("v2.2.2"),
            Version("v2.1.4"),
            Version("v2.0.2"),
           )

reallyold = (
             Version("graviton_gfort", tarballname="graviton.gfort.tar.gz"),
             Version("graviton_mod", tarballname="graviton_mod.tar.gz"),
            )

#things that get stored in the tarball if present
#graviton and graviton_mod_off are in the really old versions
storefolders = ["JHUGenerator", "JHUGenMELA", "AnalyticMELA", "graviton", "graviton_mod_off"]
#plus the manual which is treated separately

#files in this directory that should not be uploaded to the spin page
dontupload = [
              "uploadwebpage[.]py$",
              "[.]gitignore$",
              "[.]git$",
             ]
MCFMprecompiledcommit = "2d21a124764751c589ffcfc1ace2326ece6ebe6a"
websitebase = "http://spin.pha.jhu.edu"
#end of inputs
########################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--dry-run", help="create all the tarballs but don't actually upload", action="store_true", dest="dryrun")
    parser.add_argument("--dont-overwrite-mcfm", help="don't re-upload MCFM libraries that already exist", action="store_true")
    parser.add_argument("--dont-overwrite-tarballs", help="don't re-upload JHUGen tarballs that already exist", action="store_true")
    parser.add_argument("--edit-cadaverrc", help="if you want to edit the .cadaverrc before actually uploading", action="store_true")
    args = parser.parse_args()
    uploadwebpage(dontupload=dontupload, **args.__dict__)
