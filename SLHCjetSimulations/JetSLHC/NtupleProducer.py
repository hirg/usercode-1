import FWCore.ParameterSet.Config as cms

def NewNamedModule( process , Name , Object ):
  setattr( process , Name , Object )
  return getattr( process , Name )

process = cms.Process("JetNtupleProducer")
### event selection ###
# good vertices
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 15 && position.Rho <= 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)
process.load('L1TriggerConfig.GctConfigProducers.l1GctConfig_cfi')
# track quality filter
process.noscraping = cms.EDFilter("FilterOutScraping",
  applyfilter = cms.untracked.bool(True),
  debugOn = cms.untracked.bool(False),
  numtrack = cms.untracked.uint32(10),
  thresh = cms.untracked.double(0.25)
)

#process.Timing =cms.Service("Timing")

# HCAL noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.hbhefilter = cms.Path(process.HBHENoiseFilter)

### emulator ###
process.load('L1TriggerConfig.GctConfigProducers.l1GctConfig_cfi')

#Robyn added:
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")
process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/SimL1Emulator_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.l1extraParticles.centralBxOnly = cms.bool(False)

process.load("L1TriggerConfig.GctConfigProducers.l1GctConfig_cfi")
#Dump the config info
process.load("L1TriggerConfig.GctConfigProducers.l1GctConfigDump_cfi")

####
# Andy's jets
###
process.load("SLHCUpgradeSimulations.L1CaloTrigger.SLHCCaloTrigger_cff")

RingSubtractionMethods = [ "Mean" , "Median" ]# , "Constant" ]

TowerJetSizes = [ 8 , 9 , 10 , 12 ]
TowerJetShapes = [ "Circle" , "Square" ]

process.L1RingSubtractionPermutations = cms.Sequence()

for RingSubtractionMethod in RingSubtractionMethods:
  process.L1RingSubtractionPermutations += NewNamedModule(
    process , 
    RingSubtractionMethod+"RingSubtractedTower",
    cms.EDProducer( "L1RingSubtractionProducer",
      src                 = cms.InputTag("L1CaloTowerProducer"),
      RingSubtractionType = cms.string(RingSubtractionMethod) ) 
  )
 
process.L1JetAlgoPermutations = cms.Sequence()

# Produce the L1TowerJets using L1TowerJetProducer.cc for all different algo options
for TowerJetShape in TowerJetShapes:
  for TowerJetSize in TowerJetSizes:
  
    process.L1JetAlgoPermutations += NewNamedModule(
      process , 
      "TowerJet"+TowerJetShape+str(TowerJetSize)+"FromL1CaloTower" , 
      cms.EDProducer( "L1TowerJetProducer", 
        src         = cms.InputTag("L1CaloTowerProducer"),
        JetDiameter = cms.uint32( TowerJetSize ),
        JetShape    = cms.string( TowerJetShape )
      )
    )

    for RingSubtractionMethod in RingSubtractionMethods:
      process.L1JetAlgoPermutations += NewNamedModule(
        process , 
        "TowerJet"+TowerJetShape+str(TowerJetSize)+"From"+RingSubtractionMethod+"RingSubtractedTower" , 
        cms.EDProducer( "L1TowerJetProducer",
          src         = cms.InputTag( RingSubtractionMethod+"RingSubtractedTower" ),
          JetDiameter = cms.uint32( TowerJetSize ),
          JetShape    = cms.string( TowerJetShape )
        )
      )

for TowerJetShape in TowerJetShapes:
  for TowerJetSize in TowerJetSizes:
  #create different jet options

    process.L1JetAlgoPermutations += NewNamedModule(
      process,
      "TowerJet"+TowerJetShape+str(TowerJetSize)+"FromL1CaloTower",

      cms.EDProducer( "L1TowerJetProducer",
        src         = cms.InputTag("L1CaloTowerProducer"),
        JetDiameter = cms.uint32( TowerJetSize ),
        JetShape    = cms.string( TowerJetShape )
      ) 
    )

    for RingSubtractionMethod in RingSubtractionMethods:
    #create ring subtraction options
      process.L1JetAlgoPermutations += NewNamedModule(
        process,
        "TowerJet"+TowerJetShape+str(TowerJetSize)+"From"+RingSubtractionMethod+"RingSubtractedTower", 
        cms.EDProducer( "L1TowerJetProducer",
          src         = cms.InputTag(RingSubtractionMethod+"RingSubtractedTower" ),
          JetDiameter = cms.uint32  (TowerJetSize ),
          JetShape    = cms.string  (TowerJetShape)
        )
      )
      #Create 1D filtered jet options using L1TowerJetFilter1D.cc with input of L1TowerJets
      process.L1JetAlgoPermutations += NewNamedModule(
        process,
        "TowerJetFilter1D"+TowerJetShape+str(TowerJetSize)+"FromL1CaloTower",
        cms.EDProducer("L1TowerJetFilter1D",
          src                 = cms.InputTag("TowerJet"+TowerJetShape+str(TowerJetSize)+"FromL1CaloTower"), 
          NumOfOutputJets     = cms.uint32(4),
          ComparisonDirection = cms.string("eta") 
        )
      )

      #Create 2D filtered jet options using L1TowerJetFilter2D.cc with input of TowerJetFilter1D
      process.L1JetAlgoPermutations += NewNamedModule(
        process,
        "TowerJetFilter2D"+TowerJetShape+str(TowerJetSize)+"FromL1CaloTower",
        cms.EDProducer("L1TowerJetFilter2D",
          src                 = cms.InputTag("TowerJetFilter1D"+TowerJetShape+str(TowerJetSize)+"FromL1CaloTower"),
          NumOfOutputJets     = cms.uint32(12), 
          ComparisonDirection = cms.string("phi")
        )
      )

#####
# Sim ecal and hcal digis
process.load("SimCalorimetry.EcalTrigPrimProducers.ecalTrigPrimESProducer_cff")
process.load('SimCalorimetry/HcalTrigPrimProducers/hcalTTPDigis_cff')
process.simEcalTriggerPrimitiveDigis = cms.EDProducer("EcalTrigPrimProducer",
    BarrelOnly = cms.bool(False),
    InstanceEB = cms.string('ebDigis'),
    InstanceEE = cms.string('eeDigis'),
    binOfMaximum = cms.int32(6),
    Famos = cms.bool(False),
    TcpOutput = cms.bool(False),
    Debug = cms.bool(False),
    Label = cms.string('ecalDigis')
)
process.simHcalTriggerPrimitiveDigis = cms.EDProducer("HcalTrigPrimDigiProducer",
    latency = cms.int32(1),
    weights = cms.vdouble(1.0, 1.0), ##hardware algo
    #vdouble weights = { -1, -1, 1, 1} //low lumi algo
    peakFilter = cms.bool(True),
    # Input digi label (_must_ be without zero-suppression!)
    numberOfSamples = cms.int32(4),
    numberOfPresamples = cms.int32(2),
    inputLabel = cms.VInputTag(cms.InputTag('hcalDigis'),cms.InputTag('hcalDigis')),
    FG_threshold = cms.uint32(12), ## threshold for setting fine grain bit
    ZS_threshold = cms.uint32(1), ## threshold for setting fine grain bit
    MinSignalThreshold = cms.uint32(0), # For HF PMT veto
    PMTNoiseThreshold = cms.uint32(0), # For HF PMT veto
    RunZS = cms.bool(False),
    FrontEndFormatError = cms.bool(False), # Front End Format Error, for real data only
#added:
  InputTagFEDRaw = cms.InputTag("source"),
    PeakFinderAlgorithm = cms.int32(2)

)

####
# PF jet rho

process.ak5PFJets = cms.EDProducer("FastjetJetProducer",
     Active_Area_Repeats = cms.int32(1),
    doAreaFastjet = cms.bool(True),
    Ghost_EtaMax = cms.double(5.0),
    doAreaDiskApprox = cms.bool(False),
    jetType = cms.string('PFJet'),
    minSeed = cms.uint32(14327),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(True),#change from false
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(True),
    doOutputJets = cms.bool(True),
    src = cms.InputTag("towerMaker"),
    inputEtMin = cms.double(0.3),
    puPtMin = cms.double(10),
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.5)
)


##L1 AK5 jets
process.ProducePFcandsFromL1TriggerTowers = cms.EDProducer('L1CaloCandidate',
                                         srcEcal = cms.InputTag("simEcalTriggerPrimitiveDigis"),#Refers to process above
                                         srcHcal = cms.InputTag("simHcalTriggerPrimitiveDigis"),#Refers to process above
)

process.L1ak5PFJets= cms.EDProducer("FastjetJetProducer",
    Active_Area_Repeats = cms.int32(1),
    doOutputJets = cms.bool(True),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(False),
    minSeed = cms.uint32(14327),
    Ghost_EtaMax = cms.double(5.0),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(True),
    srcPVs = cms.InputTag(""),
    inputEtMin = cms.double(0.0),
    doAreaFastjet = cms.bool(True),
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    jetPtMin = cms.double(3.0),
    jetType = cms.string('PFJet'),
    src = cms.InputTag("ProducePFcandsFromL1TriggerTowers"),
    doPUOffsetCorr = cms.bool(False),
    radiusPU = cms.double(0.5),
    doAreaDiskApprox = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam       =  cms.double(0.5)
)


####
# Calo jet rho
process.ak5CaloJets = cms.EDProducer("FastjetJetProducer",
     Active_Area_Repeats = cms.int32(1),
    doAreaFastjet = cms.bool(True),
    Ghost_EtaMax = cms.double(5.0),
    doAreaDiskApprox = cms.bool(False),
    jetType = cms.string('CaloJet'),
    minSeed = cms.uint32(14327),
    voronoiRfact = cms.double(-0.9),
    doRhoFastjet = cms.bool(True),#change from false
    nSigmaPU = cms.double(1.0),
    GhostArea = cms.double(0.01),
    Rho_EtaMax = cms.double(4.4),
    useDeterministicSeed = cms.bool(True),
    doPVCorrection = cms.bool(True),
    doOutputJets = cms.bool(True),
    src = cms.InputTag("towerMaker"),
    inputEtMin = cms.double(0.3),
    puPtMin = cms.double(10),
    srcPVs = cms.InputTag("offlinePrimaryVertices"),
    jetPtMin = cms.double(3.0),
    radiusPU = cms.double(0.5),
    doPUOffsetCorr = cms.bool(False),
    inputEMin = cms.double(0.0),
    maxRecoveredHcalCells = cms.uint32(9999999),
    maxRecoveredEcalCells = cms.uint32(9999999),
    maxProblematicEcalCells = cms.uint32(9999999),
    maxBadHcalCells = cms.uint32(9999999),
    maxBadEcalCells = cms.uint32(9999999),
    maxProblematicHcalCells = cms.uint32(9999999),
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.5)
)


#Refer to these in the analyzer which is going to make the ntuple
process.demo = cms.EDAnalyzer('NtupleProducer',
  #l1 jets from l1extra particles
  extra = cms.VInputTag("l1extraParticles:Central", "l1extraParticles:Forward", "l1extraParticles:Tau"),
  #offline calo jets
  CaloJets = cms.InputTag("ak5CaloJets"),
  PFJets = cms.InputTag("ak5PFJets"),
  #Jet corrector for calo and PF jets
  JetCorrectorCalo = cms.string('ak5CaloL1FastL2L3Residual'), #jet corrections residual for data
  JetCorrectorPF = cms.string('ak5PFL1FastL2L3Residual'), #jet corrections residual for data
  RhoCaloJets = cms.InputTag("ak5CaloJets","rho"),
  RhoPFJets = cms.InputTag("ak5PFJets","rho"),
  #jet algorithms: None
  Circle8=cms.InputTag('TowerJetCircle8FromL1CaloTower'),
  Circle9=cms.InputTag('TowerJetCircle9FromL1CaloTower'),
  Circle10=cms.InputTag('TowerJetCircle10FromL1CaloTower'),
  Circle12=cms.InputTag('TowerJetCircle12FromL1CaloTower'),
  Square8=cms.InputTag('TowerJetSquare8FromL1CaloTower'),
  Square9=cms.InputTag('TowerJetSquare9FromL1CaloTower'),
  Square10=cms.InputTag('TowerJetSquare10FromL1CaloTower'),
  Square12=cms.InputTag('TowerJetSquare12FromL1CaloTower'),
  #jet algorithms: Median
  Circle8Median=cms.InputTag('TowerJetCircle8FromMedianRingSubtractedTower'),
  Circle9Median=cms.InputTag('TowerJetCircle9FromMedianRingSubtractedTower'),
  Circle10Median=cms.InputTag('TowerJetCircle10FromMedianRingSubtractedTower'),
  Circle12Median=cms.InputTag('TowerJetCircle12FromMedianRingSubtractedTower'),
  Square8Median=cms.InputTag('TowerJetSquare8FromMedianRingSubtractedTower'),
  Square9Median=cms.InputTag('TowerJetSquare9FromMedianRingSubtractedTower'),
  Square10Median=cms.InputTag('TowerJetSquare10FromMedianRingSubtractedTower'),
  Square12Median=cms.InputTag('TowerJetSquare12FromMedianRingSubtractedTower'),
  #jet algorithms: Mean
  Circle8Mean=cms.InputTag('TowerJetCircle8FromMeanRingSubtractedTower'),
  Circle9Mean=cms.InputTag('TowerJetCircle9FromMeanRingSubtractedTower'),
  Circle10Mean=cms.InputTag('TowerJetCircle10FromMeanRingSubtractedTower'),
  Circle12Mean=cms.InputTag('TowerJetCircle12FromMeanRingSubtractedTower'),
  Square8Mean=cms.InputTag('TowerJetSquare8FromMeanRingSubtractedTower'),
  Square9Mean=cms.InputTag('TowerJetSquare9FromMeanRingSubtractedTower'),
  Square10Mean=cms.InputTag('TowerJetSquare10FromMeanRingSubtractedTower'),
  Square12Mean=cms.InputTag('TowerJetSquare12FromMeanRingSubtractedTower'),
  #L1 AK5 jets
  ReRecoPfJets=cms.InputTag('L1ak5PFJets'),
  RhoL1PFJets=cms.InputTag('L1ak5PFJets:rho'),
  #Filtered jets
  FilteredCircle8=cms.InputTag('TowerJetFilter2DCircle8FromL1CaloTower'),
  FilteredCircle9=cms.InputTag('TowerJetFilter2DCircle9FromL1CaloTower'),
  FilteredCircle10=cms.InputTag('TowerJetFilter2DCircle10FromL1CaloTower'),
  FilteredCircle12=cms.InputTag('TowerJetFilter2DCircle12FromL1CaloTower'),
  FilteredSquare8=cms.InputTag('TowerJetFilter2DSquare8FromL1CaloTower'),
  FilteredSquare9=cms.InputTag('TowerJetFilter2DSquare9FromL1CaloTower'),
  FilteredSquare10=cms.InputTag('TowerJetFilter2DSquare10FromL1CaloTower'),
  FilteredSquare12=cms.InputTag('TowerJetFilter2DSquare12FromL1CaloTower'),

  fileName= cms.string('Ntuple.root'),
)


### redefine path with the pieces we want ###
process.p = cms.Path(     
    process.simEcalTriggerPrimitiveDigis
    +process.simHcalTriggerPrimitiveDigis
    +process.L1CaloTowerProducer
    +process.L1RingSubtractionPermutations
    +process.L1JetAlgoPermutations
    +process.ak5CaloJets
    +process.ak5PFJets
    +process.ProducePFcandsFromL1TriggerTowers
    +process.L1ak5PFJets
    +process.demo
)

### reconstruction
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)

# jet corrections
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.p
)


# global tag
process.GlobalTag.globaltag = 'GR_R_44_V14::All'

# N events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.Timing =cms.Service("Timing")

# input
process.source = cms.Source ("PoolSource",

   # RECO files
   fileNames = cms.untracked.vstring(
     #'file:testInput.root'
     'dcap://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/data/Run2011B/L1JetHPF/RECO/PromptReco-v1/000/179/828/86AEB242-1902-E111-B17C-003048D2C0F0.root'
     ),
   
   # RAW files
   secondaryFileNames = cms.untracked.vstring(
     'dcap://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/data/Run2011B/L1JetHPF/RAW/v1/000/179/828/589DB00F-2CFF-E011-9CCB-BCAEC5329727.root'
   ) 
)

process.TFileService=cms.Service("TFileService",
    fileName=cms.string('JetSimulations.root')
)
#


