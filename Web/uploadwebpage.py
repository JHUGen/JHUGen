#Inputs are down at the bottom.  No need to change anything above that unless there's a bug


from contextlib import contextmanager
from getpass import getpass
import os
from subprocess import check_call
import tempfile

gitdir = None

Webdir = os.getcwd()
WebGeneratordir = os.path.join(Webdir, "Generator")
if not os.path.isdir(WebGeneratordir):
     raise OSerror("You should run this from the Web directory")

def uploadwebpage():
    JHED = raw_input("Enter your JHED ID: ")
    password = getpass("Enter your JHED password: ")
    with cd("Generator"):
        create_Download(*versions)
    repmap = {
              "JHED": JHED,
              "password": password,
              "machine": "https://webdav.johnshopkins.edu/",
              "uploadfiles": getuploadfiles("."),
             }
    with NetRC(**repmap), CadaverRC(**repmap):
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

def getuploadfiles(dir):
    files = os.listdir(dir)
    result = ""
    for f in files:
        if f in dontupload:
            continue
        if os.path.isfile(os.path.join(dir, f)):
            result += "put "+os.path.join(dir, f)+"\n"
        if os.path.isdir(os.path.join(dir, f)):
            result += "mkdir {}\ncd {}\n{}cd ..\n".format(f, f, getuploadfiles(os.path.join(dir, f)))
    return result

class Version(object):
    def __init__(self, version, gitcommit=None, tarballname=None, manualname=None, manualcommit=None, melacommit=None):
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
    if not os.path.exists(os.path.join(Webdir, "Manual.pdf")):
        mostrecentversion.createtarball(force=True)
        check_call(["cp", os.path.join(gitdir, mostrecentversion.manualname), os.path.join(Webdir, "Manual.pdf")])
    for version in [mostrecentversion] + list(olderversions) + list(reallyold):
        version.createtarball()
    Download = Download_template.format(
                                        latest = mostrecentversion.getlink(),
                                        older = "\n    ".join(v.getlink() for v in [mostrecentversion] + list(olderversions)),
                                        reallyold = "\n    ".join(v.getlink() for v in reallyold),
                                       )
    with open(os.path.join(WebGeneratordir, "Download.html"), "w") as f:
        f.write(Download)

    #download MCFM libraries
    libmcfmfilename = os.path.join(WebGeneratordir, "libmcfm_7p0.so")
    if not os.path.exists(libmcfmfilename):
        check_call(["wget", "-O", libmcfmfilename, "https://github.com/JHUGen/JHUGen-backup/raw/081536e660a6712d73721eb2f7dfaded1333525f/JHUGenMELA/ggZZ_MCFM/libmcfm_7p0.so"])

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
    Compiled MCFM libraries to interface with MELA that work for a good variety of software releases: <a href="libmcfm_7p0.so">libmcfm_7p0.so</a><br>
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
versions = (
            Version("v7.0.0", gitcommit="v7.0.0.beta1", manualcommit="388d95a", melacommit="v6.9.8"),
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
              "uploadwebpage.py",
              ".gitignore",
             ]
#end of inputs
########################################################################################

if __name__ == "__main__":
    uploadwebpage()
