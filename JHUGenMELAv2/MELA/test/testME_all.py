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

referencefiles = os.listdir("reference")

for j, ref in enumerate(referencefiles):
    match = re.match("(testME_.*_Ping)_?(.*).ref", ref)
    if not match:
        raise RuntimeError("Unknown filename {}".format(ref))
    functionname = match.group(1)
    arguments = match.group(2)
    arguments = arguments.replace("2l2l", "0")
    arguments = arguments.replace("4l", "1")
    if arguments:
        arguments = [int(_) for _ in arguments.split("_")]
    else:
        arguments = []

    if functionname == "testME_ProdDec_MCFM_JHUGen_WBFZZ_Comparison_Ping":
        functionname = "testME_ProdDec_MCFM_JHUGen_WBFZZWW_Comparison_Ping"
        arguments.insert(3, 1)
    if functionname == "testME_ProdDec_MCFM_JHUGen_WBFWW_Comparison_Ping":
        functionname = "testME_ProdDec_MCFM_JHUGen_WBFZZWW_Comparison_Ping"
        arguments.insert(3, 2)

    if functionname == "testME_ProdDec_MCFM_JHUGen_JJQCDZZ_Comparison_Ping":
        functionname = "testME_ProdDec_MCFM_JHUGen_JJQCDZZWW_Comparison_Ping"
        arguments.insert(1, 1)
    if functionname == "testME_ProdDec_MCFM_JHUGen_JJQCDWW_Comparison_Ping":
        functionname = "testME_ProdDec_MCFM_JHUGen_JJQCDZZWW_Comparison_Ping"
        arguments.insert(1, 2)

    function = getattr(ROOT, functionname)
    function(*arguments)

    newfile = ref.replace(".ref", ".out")

    with open(newfile) as f:
        content = f.read()
        content = re.sub("0x[0-9a-f]*", "0xPOINTER", content)
    with open(newfile, 'w') as f:
        f.write(content)

    with open("reference/"+ref) as f:
        refcontent = f.read()
        refcontent = re.sub("0x[0-9a-f]*", "0xPOINTER", refcontent)

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
            if ratio != 0 and ratio != 1:
                badbkgratios.add(newfile)
        assert i > 1

    print j, "/", len(referencefiles)

errors = []
if different:
    errors.append("The following files are not the same as the reference files:\n    "+"\n    ".join(sorted(different)))
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
    errors.append("The following files have bad ratios between manual and automatic background:\n    "+"\n    ".join(sorted(badbkgratios)))

if errors:
    raise RuntimeError("\n\n\n".join(errors))
