import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


def addSVFit(process, cuts=None, outTableName='SVFit', path=None, USEPAIRMET=False, COMPUTEUPDOWNSVFIT=False, IsMC=False):
    srcMETTag = None
    COMPUTEMETUPDOWNSVFIT = COMPUTEUPDOWNSVFIT
    if USEPAIRMET:
        srcMETTag = cms.InputTag("corrMVAMET") if (IsMC and APPLYMETCORR) else cms.InputTag("MVAMET", "MVAMET")
    else:
        # MET corrected for central TES and EES shifts of the taus
        srcMETTag = cms.InputTag("ShiftMETcentral")

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
                                           computeForUpDownTES = cms.bool(COMPUTEUPDOWNSVFIT),
                                           computeForUpDownMET = cms.bool(COMPUTEMETUPDOWNSVFIT),
                                           METdxUP    = cms.InputTag("ShiftMETforTES", "METdxUP"),
                                           METdyUP    = cms.InputTag("ShiftMETforTES", "METdyUP"),
                                           METdxDOWN  = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
                                           METdyDOWN  = cms.InputTag("ShiftMETforTES", "METdyDOWN"),
                                           METdxUP_EES   = cms.InputTag("ShiftMETforEES", "METdxUPEES"),
                                           METdyUP_EES   = cms.InputTag("ShiftMETforEES", "METdyUPEES"),
                                           METdxDOWN_EES = cms.InputTag("ShiftMETforEES", "METdxDOWNEES"),
                                           METdyDOWN_EES = cms.InputTag("ShiftMETforEES", "METdyDOWNEES")
    )

    process.SVFitTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("SVllCandTable"),
        cut = cms.string(""),
        name = cms.string(outTableName),
        singleton = cms.bool(True), # the number of entries is variable
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
    process.svfitTask = cms.Task(process.SVllCandTable,
                                 process.SVFitTable)

    if path is None:
        process.schedule.associate(process.svfitTask)
    else:
        getattr(process, path).associate(process.svfitTask)
