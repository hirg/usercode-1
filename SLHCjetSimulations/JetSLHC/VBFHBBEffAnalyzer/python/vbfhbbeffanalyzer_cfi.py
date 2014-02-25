import FWCore.ParameterSet.Config as cms

vbfHBBEffAnalyzer = cms.EDAnalyzer('VBFHBBEffAnalyzer',
                                   JetsFromHiggsDecay=cms.bool(True),
                                   JetsFromVBF=cms.bool(False),
                                   FilteredCircle8= cms.InputTag("L1TowerJetFilter2D"),
                                   extra = cms.VInputTag("l1extraParticles:Central", "l1extraParticles:Forward", "l1extraParticles:Tau"),
                                   UWJets= cms.InputTag("L1CaloJetExpander")
)
