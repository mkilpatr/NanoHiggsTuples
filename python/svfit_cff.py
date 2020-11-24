import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *


def addSVFit(process, cuts=None, outTableName='SVFit', path=None, USEPAIRMET=False, COMPUTEUPDOWNSVFIT=False):
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

    process.svfitTask = cms.Task(
        process.SVllCandTable,
        )

    if path is None:
        process.schedule.associate(process.svfitTask)
    else:
        getattr(process, path).associate(process.svfitTask)
