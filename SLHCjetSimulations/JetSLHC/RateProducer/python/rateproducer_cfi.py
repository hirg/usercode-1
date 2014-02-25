import FWCore.ParameterSet.Config as cms

rateProducer = cms.EDAnalyzer('RateProducer',
                                   FilteredCircle8= cms.InputTag("L1TowerJetFilter2D"),
                                   extra = cms.VInputTag("l1extraParticles:Central", "l1extraParticles:Forward", "l1extraParticles:Tau"),
                                   UWJets= cms.InputTag("L1CaloJetExpander")
)
