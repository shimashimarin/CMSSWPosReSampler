# coding: utf-8

import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# setup minimal options
#options = VarParsing("python")
#options.setDefault("inputFiles", "root://xrootd-cms.infn.it//store/mc/RunIIFall17MiniAOD/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v2/00000/9A439935-1FFF-E711-AE07-D4AE5269F5FF.root")  # noqa
#options.parseArguments()

# define the process to run
process = cms.Process("Uwt")

# minimal configuration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.source = cms.Source("PoolSource",
   fileNames=cms.untracked.vstring('file:/home/pku/stqian/NegWgt/CMSSW_11_3_0_pre2/src/RwtUwt/RwtUwt/python/myOutputFile.root'))
# process.source = cms.Source("EmptySource")
# process options
process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
)

process.UwtFilter = cms.EDFilter("UwtFilter", ifkeep = cms.InputTag("RwtUwt","ifkeep"))

# setup MyPlugin by loading the auto-generated cfi (see MyPlugin.fillDescriptions)
#process.load("XGB_Example.XGBoostExample.XGBoostExample_cfi")
# process.XGBoostExample.model_path = cms.string("/Your/Path/data/lowVer.model")
# process.XGBoostExample.test_data_path = cms.string("/Your/Path/data/Test_data.csv")

# define what to run in the path

process.out = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    ),
    fileName = cms.untracked.string('test_filtering.root')
)

process.maxEvents.input = 2000

process.p1 = cms.Path(process.UwtFilter)
process.e = cms.EndPath(process.out)
