import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

# Load tools and function definitions
from PhysicsTools.SmuontorUtils.tools.vid_id_tools import *

def addMuon(process, cuts=None, outTableName='Muon', path=None, IsMC=False):
    ### Mu Ghost cleaning
    process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                       src = cms.InputTag("slimmedMuons"),
                                       presmuontion = cms.string("track.isNonnull"),
                                       passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                       fractionOfSharedSegments = cms.double(0.499))
    

    process.bareSoftMuons = cms.EDFilter("PATMuonRefSmuontor",
        src = cms.InputTag("cleanedMu"),
        cut = cms.string("isLooseMuon && pt>10")
    )
    
    process.softMuons = cms.EDProducer("MuFiller",
        src = cms.InputTag("bareSoftMuons"),
        genCollection = cms.InputTag("prunedGenParticles"),
        rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
        vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
        sampleType = cms.int32(LEPTON_SETUP),                     
        setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
        cut = cms.string(""),
        flags = cms.PSet(
            ID = cms.string("userFloat('isPFMuon')" ), # PF ID
            isGood = cms.string(MUCUT)
        )
    )
    
    process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons+ process.softMuons)
    
    process.muonTask = cms.Task(process.muons)

    if path is None:
        process.schedule.associate(process.muonTask)
    else:
        getattr(process, path).associate(process.muonTask)
