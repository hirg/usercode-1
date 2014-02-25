import FWCore.ParameterSet.Config as cms

calibJetProducer = cms.EDProducer('CalibTowerJetProducer',
                                  FilteredCircle8=cms.InputTag('TowerJetFilter2DCircle8FromL1CaloTower'),  
)
