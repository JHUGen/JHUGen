#!/usr/bin/env python

import os
import re

import ROOT

ROOT.gROOT.Macro("loadMELA.C")
ROOT.gROOT.LoadMacro("testME_v2.c+")

for ref in os.listdir("reference"):
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
