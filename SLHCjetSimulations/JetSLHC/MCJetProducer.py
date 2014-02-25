import FWCore.ParameterSet.Config as cms

def NewNamedModule( process , Name , Object ):
  setattr( process , Name , Object )
  return getattr( process , Name )

process = cms.Process("UUU")
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

process.Timing =cms.Service("Timing")

# HCAL noise filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.hbhefilter = cms.Path(process.HBHENoiseFilter)

### emulator ###
process.load('L1TriggerConfig.GctConfigProducers.l1GctConfig_cfi')

#Robyn added:
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cfi")
process.load("Configuration.StandardSequences.RawToDigi_cff")
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


### redefine path with the pieces we want ###
process.p = cms.Path(     
    process.simEcalTriggerPrimitiveDigis
    +process.simHcalTriggerPrimitiveDigis
    +process.L1CaloTowerProducer
    +process.L1RingSubtractionPermutations
    +process.L1JetAlgoPermutations
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
     'file:/data/pioppi/work/SLHC/SIGNAL/src/0485535E-4001-E111-9D12-E0CB4E19F96C.root'
   ),
   
)

process.TFileService=cms.Service("TFileService",
    fileName=cms.string('MCJetSimulations.root')
)
