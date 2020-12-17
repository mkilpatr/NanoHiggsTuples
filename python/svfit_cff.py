import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

def addSVFit(process, cuts=None, outTableName='SVFit', path=None, USEPAIRMET=False, COMPUTEUPDOWNSVFIT=False, IsMC=False):
    srcMETTag = None
    if USEPAIRMET:
        srcMETTag = cms.InputTag("corrMVAMET") if (IsMC and APPLYMETCORR) else cms.InputTag("MVAMET", "MVAMET")
    else:
        # MET corrected for central TES and EES shifts of the taus
        srcMETTag = cms.InputTag("ShiftMETcentral")


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
				 process.SVllCandTable,
                                 process.SVFitTable)

    if path is None:
        process.schedule.associate(process.svfitTask)
    else:
        getattr(process, path).associate(process.svfitTask)
