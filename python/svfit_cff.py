import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def addMETProcesses(process, cuts=None, path=None, USEPAIRMET=True, COMPUTEUPDOWNSVFIT=False, IsMC=False, APPLYMETCORR=True):
    process.METSequence = cms.Sequence()
    if USEPAIRMET:
        print "Using pair MET (MVA MET)"
        from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
        runMVAMET(process, jetCollectionPF = "patJetsReapplyJEC")
        process.MVAMET.srcLeptons = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
        process.MVAMET.requireOS = cms.bool(False)
        process.MVAMET.permuteLeptonsWithinPlugin = cms.bool(False)
        process.MVAMET.leptonPermutations = cms.InputTag("barellCand")
    
        process.MVAMETInputs = cms.Sequence(
            process.slimmedElectronsTight + process.slimmedMuonsTight + process.slimmedTausLoose + process.slimmedTausLooseCleaned + process.patJetsReapplyJECCleaned +
            process.pfCHS + process.pfChargedPV + process.pfChargedPU + process.pfNeutrals + process.neutralInJets +
            process.pfMETCands + process.pfTrackMETCands + process.pfNoPUMETCands + process.pfPUCorrectedMETCands + process.pfPUMETCands +
            process.pfChargedPUMETCands + process.pfNeutralPUMETCands + process.pfNeutralPVMETCands + process.pfNeutralUnclusteredMETCands +
            process.pfChs +
            process.ak4PFCHSL1FastjetCorrector + process.ak4PFCHSL2RelativeCorrector + process.ak4PFCHSL3AbsoluteCorrector + process.ak4PFCHSResidualCorrector +
            process.ak4PFCHSL1FastL2L3Corrector + process.ak4PFCHSL1FastL2L3ResidualCorrector +
            process.tauDecayProducts + process.tauPFMET + process.tauMET + process.tausSignificance
        )
        for met in ["pfMET", "pfTrackMET", "pfNoPUMET", "pfPUCorrectedMET", "pfPUMET", "pfChargedPUMET", "pfNeutralPUMET", "pfNeutralPVMET", "pfNeutralUnclusteredMET"]:
            process.MVAMETInputs += getattr(process, met)
            process.MVAMETInputs += getattr(process, "ak4JetsFor"+met)
            process.MVAMETInputs += getattr(process, "corr"+met)
            process.MVAMETInputs += getattr(process, met+"T1")
            process.MVAMETInputs += getattr(process, "pat"+met)
            process.MVAMETInputs += getattr(process, "pat"+met+"T1")        
    
        process.METSequence += cms.Sequence(process.MVAMETInputs + process.MVAMET)
    
    
    else:
        print "Using event pfMET (same MET for all pairs)"
    
        PFMetName = "slimmedMETs"
        uncorrPFMetTag = cms.InputTag(PFMetName)
    
        # patch to get a standalone MET significance collection
        process.METSignificance = cms.EDProducer ("ExtractMETSignificance",
                                                      #srcMET=cms.InputTag(PFMetName,"","TEST")
                                                      srcMET=uncorrPFMetTag
                                                      )
    
        # Shift met due to central corrections of TES and EES
        process.ShiftMETcentral = cms.EDProducer ("ShiftMETcentral",
                                                  srcMET = uncorrPFMetTag,
                                                  tauUncorrected = cms.InputTag("bareTaus"),
                                                  tauCorrected = cms.InputTag("softTaus")
                                                  )
    
    
        process.METSequence += process.METSignificance
        process.METSequence += process.ShiftMETcentral

    ## ----------------------------------------------------------------------
    ## Z-recoil correction
    ## ----------------------------------------------------------------------
    
    # corrMVAPairMET = []
    if IsMC and APPLYMETCORR:
        if USEPAIRMET:
            process.selJetsForZrecoilCorrection = cms.EDFilter("PATJetSelector",
                src = cms.InputTag("jets"),                                      
                cut = cms.string("pt > 30. & abs(eta) < 4.7"), 
                filter = cms.bool(False)
            )
            process.corrMVAMET = cms.EDProducer("ZrecoilCorrectionProducer",                                                   
                srcPairs = cms.InputTag("barellCand"),
                srcMEt = cms.InputTag("MVAMET", "MVAMET"),
                srcGenParticles = cms.InputTag("prunedGenParticles"),
                srcJets = cms.InputTag("selJetsForZrecoilCorrection"),
                correction = cms.string("HTT-utilities/RecoilCorrections/data/MvaMET_MG_2016BCD.root")
            )
            process.METSequence += process.selJetsForZrecoilCorrection        
            process.METSequence += process.corrMVAMET
    
        else:
            raise ValueError("Z-recoil corrections for PFMET not implemented yet !!")
    
    
    srcMETTag = None
    if USEPAIRMET:
      srcMETTag = cms.InputTag("corrMVAMET") if (IsMC and APPLYMETCORR) else cms.InputTag("MVAMET", "MVAMET")
    else:
      # MET corrected for central TES and EES shifts of the taus
      srcMETTag = cms.InputTag("ShiftMETcentral")

    return process

def addSVFit(process, cuts=None, outTableName='SVFit', path=None, USEPAIRMET=True, COMPUTEUPDOWNSVFIT=False, IsMC=False, APPLYMETCORR=True):
    addMETProcesses(process, USEPAIRMET=USEPAIRMET, COMPUTEUPDOWNSVFIT=COMPUTEUPDOWNSVFIT, IsMC=IsMC, APPLYMETCORR=APPLYMETCORR)
    #Leptons
    muString = "softMuons"
    eleString = "slimmedElectrons"#"softElectrons"
    tauString = "softTaus"
    process.softLeptons = cms.EDProducer("CandViewMerger",
        #src = cms.VInputTag(cms.InputTag("slimmedMuons"), cms.InputTag("slimmedElectrons"),cms.InputTag("slimmedTaus"))
        src = cms.VInputTag(cms.InputTag(muString), cms.InputTag(eleString),cms.InputTag(tauString))
    )

    BUILDONLYOS = True
    LLCUT="mass>0"
    ##
    ## Build ll candidates (here OS)
    ##
    decayString="softLeptons softLeptons"
    checkcharge=False
    if BUILDONLYOS:
        decayString="softLeptons@+ softLeptons@-"
        checkcharge=True
    process.barellCand = cms.EDProducer("CandViewShallowCloneCombiner",
                                        decay = cms.string(decayString),
                                        cut = cms.string(LLCUT),
                                        checkCharge = cms.bool(checkcharge)
    )
    
    ## ----------------------------------------------------------------------
    ## SV fit
    ## ----------------------------------------------------------------------
    #if USECLASSICSVFIT:
    #    print "Using CLASSIC_SV_FIT"
    process.SVllCandTable = cms.EDProducer("ClassicSVfitInterface",
                                           srcPairs   = cms.InputTag("barellCand"),
                                           srcSig     = cms.InputTag("METSignificance", "METSignificance"),
                                           srcCov     = cms.InputTag("METSignificance", "METCovariance"),
                                           usePairMET = cms.bool(USEPAIRMET),
                                           srcMET     = srcMETTag,
					   SVFitName  = cms.string(outTableName),
                                           debug      = cms.bool(True),
    )

    process.SVFitTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("SVllCandTable"),
        cut = cms.string(""),
        name = cms.string(outTableName),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(True),
        variables = cms.PSet( P4Vars,
            dz    = Var("dz()", float, doc = "pfcand info dz", precision=8),
            fromPV= Var("fromPV()", float, doc = "pfcand info from Primary Vertex", precision=8),
        ),
    )

    process.SVFitTable.variables.pt.precision=10
    process.SVFitTable.variables.eta.precision=12
    process.SVFitTable.variables.phi.precision=10
    process.SVFitTable.variables.mass.precision=10
    process.svfitTask = cms.Task(process.softLeptons,
				 process.barellCand, 
				 process.ShiftMETcentral,
				 process.SVllCandTable,
                                 process.SVFitTable)

    if path is None:
        process.schedule.associate(process.svfitTask)
    else:
        getattr(process, path).associate(process.svfitTask)
