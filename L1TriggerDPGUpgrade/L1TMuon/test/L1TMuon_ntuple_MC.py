##################################
#  L1TMuon_ntuple_MC.py
#  - python module that runs /plugins/rpc_tuplizer.cc
#  - creates ntuple of L1TMuon trigger primitives and tracks
#  - currently built for Montel Carlo
#
#  - written by David Curry
##################################

import FWCore.ParameterSet.Config as cms

# Load all process
process = cms.Process('L1TMUONtuple')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1TMuonTriggerPrimitiveProducer_cfi')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1CSCTFTrackConverter_cfi')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1DTTFTrackConverter_cfi')
process.load('L1TriggerDPGUpgrade.L1TMuon.L1RPCTFTrackConverter_cfi')

# gen particles - its one of these
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


# Global Tag to specify the "good" events to run over
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

# Import Monte Carlo data file
infile = ['file:/exports/uftrig01a/dcurry/data/rpc/SingleMuLowPt_5GeVto200GeV_GEN_SIM_DIGI_L1_1_20_2015.root']

#infile = ['file:/exports/uftrig01a/dcurry/data/monte_carlo/neutrinoGun_2012_GEN_SIM_DIGI_RAW_003048F0E186.root']

process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring(infile) )

# How many events to run over
process.maxEvents = cms.untracked.PSet( input=cms.untracked.int32(100) )

# Display event progress during runtime
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# ------- Main Process: rpc_tuplizer.cc -----------

process.tuple = cms.EDAnalyzer("rpc_tuplizer",
                               genSrc       = cms.untracked.InputTag("genParticles"),
                               muonsTag     = cms.InputTag("muons"),
                               csctfTag     = cms.InputTag("L1CSCTFTrackConverter"),
                               rpcTPTag     = cms.InputTag("L1TMuonTriggerPrimitives","RPC"),
                               cscSegTag    = cms.InputTag("cscSegments"),
                               outputFile   = cms.string("/exports/uftrig01a/dcurry/data/rpc/singleMu_lowPt_clustering_3_3_15.root"),
                               printLevel   = cms.untracked.int32(2),
                               isMC         = cms.untracked.int32(1),
                               lutParam     = cms.PSet(
                                                    isBeamStartConf = cms.untracked.bool(True),
                                                    ReadPtLUT = cms.bool(False),
                                                    PtMethod = cms.untracked.uint32(32)
                                                    )
                               )


#--------------------------------------------------

# ------- L1TMuon Modules and final process -------------------------

process.L1TMuonSeq = cms.Sequence( process.genParticles             +
                                   process.L1TMuonTriggerPrimitives +
                                   process.L1CSCTFTrackConverter    +
                                   process.L1DTTFTrackConverter     +
                                   process.L1RPCTFTrackConverters   +
                                   process.tuple
                                   )
#--------------------------------------------------

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)
process.Schedule = cms.Schedule(process.L1TMuonPath)


# ------ Uncomment below to create root file with all objects created during runtime
#process.myOut = cms.OutputModule("PoolOutputModule",
#                                 fileName = cms.untracked.string('saveEverything.root'),
#                                 outputCommands = cms.untracked.vstring('keep *') )
#process.outpath = cms.EndPath(process.myOut)




