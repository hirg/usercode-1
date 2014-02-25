import FWCore.ParameterSet.Config as cms

calibTowerJetProducer = cms.EDProducer('CalibTowerJetProducer',
                                   FilteredCircle8= cms.InputTag("L1TowerJetFilter2D"),
                                   FilteredFwdCircle8=cms.InputTag("L1TowerFwdJetFilter2D"),
                                   inPtdata_file = cms.FileInPath(
                                       "ProduceJetCollections/CalibTowerJetProducer/data/looselogval.log"),
                                   inRhodata_file = cms.FileInPath(
                                       "ProduceJetCollections/CalibTowerJetProducer/data/rho_lookup.txt"),
                                   inMVA_weights_file =  cms.FileInPath(
                                       "ProduceJetCollections/CalibTowerJetProducer/data/TMVARegression_BDT.weights.xml"),

)
