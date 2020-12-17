import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig

def addTaus(process, cuts=None, outTableName='Taus', path=None, USEPAIRMET=False, IsMC=False):
    ##
    ## Taus
    ##
    PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes

    process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
      src = cms.InputTag("offlineSlimmedPrimaryVertices"),
      cut = cms.string(PVERTEXCUT),
      filter = cms.bool(False), # if True, rejects events . if False, produce emtpy vtx collection
    )

    updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
    TAUCUT="tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 && pt>18" #miniAOD 
    APPLYTESCORRECTION=False
    YEAR="2018"
    TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"

    tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = True,
                        updatedTauName = updatedTauName,
                        toKeep = ["deepTau2017v2p1", "2017v2"]  #["2017v1", "dR0p32017v2"]
    )
    
    tauIdEmbedder.runTauID()
    
    # old sequence starts here
    process.bareTaus = cms.EDFilter("PATTauRefSelector",
       src = cms.InputTag(updatedTauName), 
       cut = cms.string(TAUCUT),
       )
    
    # TES corrections: https://indico.cern.ch/event/887196/contributions/3743090/attachments/1984772/3306737/TauPOG_TES_20200210.pdf
    
    # EES corrections: https://indico.cern.ch/event/868279/contributions/3665970/attachments/1959265/3267731/FES_9Dec_explained.pdf
    
    NominalTESCorrection=cms.double(0.0) #in percent\
    APPLYTESCORRECTION = APPLYTESCORRECTION if IsMC else False # always false if data
    
    # 2016 data - DeepTau2017v2p1
    NomTESUncDM0      = cms.double(0.8)  # in percent, up/down uncertainty of TES
    NomTESUncDM1      = cms.double(0.6)  # in percent, up/down uncertainty of TES
    NomTESUncDM10     = cms.double(0.8)  # in percent, up/down uncertainty of TES
    NomTESUncDM11     = cms.double(1.1)  # in percent, up/down uncertainty of TES
    NomTESCorDM0      = cms.double(-0.9) # DecayMode==0
    NomTESCorDM1      = cms.double(-0.1) # DecayMode==1
    NomTESCorDM10     = cms.double(0.3)  # DecayMode==10
    NomTESCorDM11     = cms.double(-0.2) # DecayMode==11
    
    #EES BARREL
    NomEFakeESCorDM0B     = cms.double(0.679) #DecayMode==0
    NomEFakeESUncDM0BUp    = cms.double(0.806) #DecayMode==0
    NomEFakeESUncDM0BDown  = cms.double(0.982) #DecayMode==0
    NomEFakeESCorDM1B      = cms.double(3.389) #DecayMode==1
    NomEFakeESUncDM1BUp    = cms.double(1.168) #DecayMode==1
    NomEFakeESUncDM1BDown  = cms.double(2.475) #DecayMode==1
    #EES ENDCAP
    NomEFakeESCorDM0E      = cms.double(-3.5)   #DecayMode==0
    NomEFakeESUncDM0EUp    = cms.double(1.808)  #DecayMode==0
    NomEFakeESUncDM0EDown  = cms.double(1.102)  #DecayMode==0
    NomEFakeESCorDM1E      = cms.double(5.)      #DecayMode==1
    NomEFakeESUncDM1EUp    = cms.double(6.57)   #DecayMode==1
    NomEFakeESUncDM1EDown  = cms.double(5.694)  #DecayMode==1
    
    TESyear = "2016Legacy"
    
    # 2017 data - DeepTau2017v2p1
    if YEAR == 2017:
        NomTESUncDM0      = cms.double(1.0)  # in percent, up/down uncertainty of TES
        NomTESUncDM1      = cms.double(0.6)  # in percent, up/down uncertainty of TES
        NomTESUncDM10     = cms.double(0.7)  # in percent, up/down uncertainty of TES
        NomTESUncDM11     = cms.double(1.4)  # in percent, up/down uncertainty of TES
        NomTESCorDM0      = cms.double(0.4)  # DecayMode==0
        NomTESCorDM1      = cms.double(0.2)  # DecayMode==1
        NomTESCorDM10     = cms.double(0.1)  # DecayMode==10
        NomTESCorDM11     = cms.double(-1.3) # DecayMode==1
    
        #EES BARREL
        NomEFakeESCorDM0B      = cms.double(0.911) #DecayMode==0
        NomEFakeESUncDM0BUp    = cms.double(1.343) #DecayMode==0
        NomEFakeESUncDM0BDown  = cms.double(0.882) #DecayMode==0
        NomEFakeESCorDM1B      = cms.double(1.154) #DecayMode==1
        NomEFakeESUncDM1BUp    = cms.double(2.162) #DecayMode==1
        NomEFakeESUncDM1BDown  = cms.double(0.973) #DecayMode==1
        #EES ENDCAP
        NomEFakeESCorDM0E      = cms.double(-2.604)   #DecayMode==0
        NomEFakeESUncDM0EUp    = cms.double(2.249)    #DecayMode==0
        NomEFakeESUncDM0EDown  = cms.double(1.43)     #DecayMode==0
        NomEFakeESCorDM1E      = cms.double(1.5)    #DecayMode==1
        NomEFakeESUncDM1EUp    = cms.double(6.461)      #DecayMode==1
        NomEFakeESUncDM1EDown  = cms.double(4.969)    #DecayMode==1
    
        TESyear = "2017ReReco"
    
    # 2018 data - DeepTau2017v2p1
    if YEAR == 2018:
        NomTESUncDM0          = cms.double(0.9)  # in percent, up/down uncertainty of TES
        NomTESUncDM1          = cms.double(0.5)  # in percent, up/down uncertainty of TES
        NomTESUncDM10         = cms.double(0.7)  # in percent, up/down uncertainty of TES
        NomTESUncDM11         = cms.double(1.2)  # in percent, up/down uncertainty of TES
        NomTESCorDM0          = cms.double(-1.6) # DecayMode==0
        NomTESCorDM1          = cms.double(-0.5) # DecayMode==1
        NomTESCorDM10         = cms.double(-1.2) # DecayMode==10
        NomTESCorDM11         = cms.double(-0.4) # DecayMode==11
    
        #EES BARREL
        NomEFakeESCorDM0B      = cms.double(1.362)    #DecayMode==0
        NomEFakeESUncDM0BUp    = cms.double(0.904)    #DecayMode==0
        NomEFakeESUncDM0BDown  = cms.double(0.474)    #DecayMode==0
        NomEFakeESCorDM1B      = cms.double(1.954)    #DecayMode==1
        NomEFakeESUncDM1BUp    = cms.double(1.226)    #DecayMode==1
        NomEFakeESUncDM1BDown  = cms.double(1.598)    #DecayMode==1
        #EES ENDCAP
        NomEFakeESCorDM0E      = cms.double(-3.097)   #DecayMode==0
        NomEFakeESUncDM0EUp    = cms.double(3.404)    #DecayMode==0
        NomEFakeESUncDM0EDown  = cms.double(1.25)     #DecayMode==0
        NomEFakeESCorDM1E      = cms.double(-1.5)     #DecayMode==1
        NomEFakeESUncDM1EUp    = cms.double(5.499)    #DecayMode==1
        NomEFakeESUncDM1EDown  = cms.double(4.309)    #DecayMode==1
    
        TESyear = "2018ReReco"
    
    process.softTaus = cms.EDProducer("TauFiller",
       src = cms.InputTag("bareTaus"),
       genCollection = cms.InputTag("prunedGenParticles"),
       vtxCollection = cms.InputTag("goodPrimaryVertices"),
       cut = cms.string(TAUCUT),
       discriminator = cms.string(TAUDISCRIMINATOR),

       NominalTESCorrection             = NominalTESCorrection,    
       NominalTESUncertaintyDM0         = NomTESUncDM0,
       NominalTESUncertaintyDM1         = NomTESUncDM1,
       NominalTESUncertaintyDM10        = NomTESUncDM10,
       NominalTESUncertaintyDM11        = NomTESUncDM11,
       NominalTESCorrectionDM0          = NomTESCorDM0,
       NominalTESCorrectionDM1          = NomTESCorDM1,
       NominalTESCorrectionDM10         = NomTESCorDM10,
       NominalTESCorrectionDM11         = NomTESCorDM11,
    
       NominalEFakeESCorrectionDM0B      = NomEFakeESCorDM0B,
       NominalEFakeESUncertaintyDM0BUp   = NomEFakeESUncDM0BUp, 
       NominalEFakeESUncertaintyDM0BDown = NomEFakeESUncDM0BDown, 
       NominalEFakeESCorrectionDM1B      = NomEFakeESCorDM1B,
       NominalEFakeESUncertaintyDM1BUp   = NomEFakeESUncDM1BUp, 
       NominalEFakeESUncertaintyDM1BDown = NomEFakeESUncDM1BDown, 
       NominalEFakeESCorrectionDM0E      = NomEFakeESCorDM0E,
       NominalEFakeESUncertaintyDM0EUp   = NomEFakeESUncDM0EUp, 
       NominalEFakeESUncertaintyDM0EDown = NomEFakeESUncDM0EDown, 
       NominalEFakeESCorrectionDM1E      = NomEFakeESCorDM1E,
       NominalEFakeESUncertaintyDM1EUp   = NomEFakeESUncDM1EUp, 
       NominalEFakeESUncertaintyDM1EDown = NomEFakeESUncDM1EDown, 
    
       ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
       # ApplyTESUpDown = cms.bool(True if IsMC else False), # no shift computation when data
       flags = cms.PSet(
            isGood = cms.string("")
            ),
    
       year = cms.string(TESyear)
       )

    process.newtaus=cms.Sequence(process.rerunMvaIsolationSequence + process.slimmedTausNewID + process.bareTaus)   
 
    process.taus=cms.Task(process.newtaus, 
                          process.softTaus)

    if path is None:
        process.schedule.associate(process.taus)
    else:
        getattr(process, path).associate(process.taus)
