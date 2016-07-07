from contextlib import contextmanager
import os
import subprocess
import tempfile

gitdir = None

class Version(object):
    store = ["JHUGenerator", "JHUGenMELA", "JHUGenMELAv2", "AnalyticMELA", "graviton", "graviton_mod_off"]
    #plus the manual which is treated separately

    def __init__(self, version, gitcommit=None, tarballname=None, manualname=None, manualcommit=None):
        self.version = version
        if gitcommit is None: gitcommit = version
        if tarballname is None: tarballname = "JHUGenerator.{}.tar.gz".format(self.version)
        if manualname is None: manualname = "manJHUGenerator.{}.pdf".format(self.version)
        if manualcommit is None: manualcommit = gitcommit
        self.gitcommit = gitcommit
        self.WebGeneratordir = os.getcwd()
        self.tarballname = tarballname
        self.manualname = manualname
        self.manualcommit = manualcommit

    def getlink(self):
        return link_template.format(tarballname=self.tarballname)

    def createtarball(self):
        global gitdir

        if gitdir is None:
            tmpdir = tempfile.mkdtemp()
            with cd(tmpdir):
                subprocess.check_call(["git", "clone", "git@github.com:hroskes/JHUGen"])
                gitdir = os.path.join(tmpdir, "JHUGen")

        with cd(gitdir):
            subprocess.check_call(["git", "clean", "-f", "-x", "-d"])
            store = []

            #manual
            subprocess.check_call(["git", "checkout", self.manualcommit])
            try:
                with cd("Manual"):
                    for i in range(3):
                        subprocess.check_call(["pdflatex.exe", "manJHUGenerator.tex"])
                    subprocess.check_call(["mv", "manJHUGenerator.pdf", "../"+self.manualname])
                    store.append(self.manualname)
            except OSError: #Manual directory doesn't exist, really old version
                pass

            subprocess.check_call(["git", "checkout", self.gitcommit])
            for folder in self.store:
                if os.path.exists(folder):
                    store.append(folder)
            subprocess.check_call(["tar", "-cvzf", self.tarballname] + store)
            subprocess.check_call(["mv", self.tarballname, self.WebGeneratordir])

def create_Download(mostrecentversion, *olderversions):
    reallyold = [
                 Version("graviton_gfort", tarballname="graviton.gfort.tar.gz"),
                 Version("graviton_mod", tarballname="graviton_mod.tar.gz"),
                ]
    mostrecentversion.createtarball()
    subprocess.check_call(["cp", os.path.join(gitdir, mostrecentversion.manualname), "../Manual.pdf"])
    for version in [mostrecentversion] + list(olderversions) + reallyold:
        version.createtarball()
    Download = Download_template.format(
                                        latest = mostrecentversion.getlink(),
                                        older = "\n    ".join(o.getlink() for o in olderversions),
                                        reallyold = "\n    ".join(r.getlink() for r in reallyold),
                                       )
    with open("Download.html", "w") as f:
        f.write(Download)

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
    Older versions:<br>
    {older}
    <br>
    Much older versions:<br>
    {reallyold}
    <br>
  </body>
</html>
""".lstrip()

