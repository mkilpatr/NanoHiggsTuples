import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

def addElec(process, cuts=None, outTableName='Elec', path=None, IsMC=False, nanosec="25"):
    #constants
    ELECUT="pt>7"#"gsfTrack.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"

    process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
    
    #**********************
    dataFormat = DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process, dataFormat)
    #**********************
    
    process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
    # overwrite a default parameter: for miniAOD, the collection name is a slimmed one
    #process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
    
    from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
    process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

    # Define which IDs we want to produce
    # Each of these two example IDs contains all four standard 
    # cut-based ID working points (only two WP of the PU20bx25 are actually used here).
    my_id_modules =[
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',    # both 25 and 50 ns cutbased ids produced
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',                 # recommended for both 50 and 25 ns
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff', # will not be produced for 50 ns, triggering still to come
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',    # 25 ns trig
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_50ns_Trig_V1_cff',    # 50 ns trig
    ] 
    #['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_'+nanosec+'ns_V1_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
    #Add them to the VID producer
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
        
    
    #process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
    #   src = cms.InputTag("slimmedElectrons"),#"calibratedPatElectrons"),
    #   cut = cms.string(ELECUT)
    #   )
    
    
    process.softElectrons = cms.EDProducer("EleFiller",
       src    = cms.InputTag("slimmedElectrons"),
       rhoCollection = cms.InputTag("fixedGridRhoFastjetAll",""),
       vtxCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
       genCollection = cms.InputTag("prunedGenParticles"),
       sampleType = cms.int32(2012),          
       setup = cms.int32(2012), # define the set of effective areas, rho corrections, etc.
    
       #CUT BASED ELE ID
       electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-veto"),
       electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-tight"),
       electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-medium"),
       electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-"+nanosec+"ns-V1-standalone-loose"),
    
       #MVA ELE ID (only for 25ns right now)
       eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
       eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
       mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
       mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
    
    #    cut = cms.string("userFloat('SIP')<100"),
    #   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),
       cut = cms.string(ELECUT),
       flags = cms.PSet(
            ID = cms.string("userInt('isBDT')"), # BDT MVA ID
            isGood = cms.string("")
            )
       )
    
    egmMod = 'egmGsfElectronIDs'
    mvaMod = 'electronMVAValueMapProducer'
    egmSeq = 'egmGsfElectronIDSequence'
    setattr(process,egmMod,process.egmGsfElectronIDs.clone())
    setattr(process,mvaMod,process.electronMVAValueMapProducer.clone())
    setattr(process,egmSeq,cms.Sequence(getattr(process,mvaMod)*getattr(process,egmMod)))
    process.electrons = cms.Sequence(getattr(process,mvaMod)*getattr(process,egmMod) * process.softElectrons)#process.bareSoftElectrons


    ### ----------------------------------------------------------------------
    ### Lepton Cleaning (clean electrons collection from muons)
    ### ----------------------------------------------------------------------
    
    process.cleanSoftElectrons = cms.EDProducer("PATElectronCleaner",
        # pat electron input source
        src = cms.InputTag("softElectrons"),
        # preselection (any string-based cut for pat::Electron)
        preselection = cms.string(''),
        # overlap checking configurables
        checkOverlaps = cms.PSet(
            muons = cms.PSet(
               src       = cms.InputTag("softMuons"), # Start from loose lepton def
               algorithm = cms.string("byDeltaR"),
               preselection        = cms.string("(isGlobalMuon || userFloat('isPFMuon'))"), #
               deltaR              = cms.double(0.05),  
               checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
               pairCut             = cms.string(""),
               requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
            )
        ),
        # finalCut (any string-based cut for pat::Electron)
        finalCut = cms.string(''),
    )


    process.elecTask = cms.Task(getattr(process,mvaMod), 
				getattr(process,egmMod), 
				process.softElectrons,
				process.cleanSoftElectrons)

    if path is None:
        process.schedule.associate(process.elecTask)
    else:
        getattr(process, path).associate(process.elecTask)
