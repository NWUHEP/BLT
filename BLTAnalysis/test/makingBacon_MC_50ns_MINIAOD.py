import FWCore.ParameterSet.Config as cms
import os

process = cms.Process('MakingBacon')

is_data_flag  = False                                      # flag for if process data
do_hlt_filter = False                                      # flag to skip events that fail relevant triggers
hlt_filename  = "BaconAna/DataFormats/data/HLTFile_50ns"   # list of relevant triggers

cmssw_base = os.environ['CMSSW_BASE']

#--------------------------------------------------------------------------------
# Import of standard configurations
#================================================================================
process.load('FWCore/MessageService/MessageLogger_cfi')
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('TrackingTools/TransientTrack/TransientTrackBuilder_cfi')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MCRUN2_74_V9A::All'

# produce photon isolation with proper footprint removal
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
#PUPPI JEC
process.load('BaconProd/Ntupler/myPUPPICorrections_cff')
process.puppijec.connect = cms.string('sqlite:////'+cmssw_base+'/src/BaconProd/Utils/data/PY8_RunIISpring15DR74_bx50_MC.db')
#CHS
process.pfNoPileUpJME = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

#--------------------------------------------------------------------------------
# Import custom configurations
#================================================================================
# custom jet stuff (incl. GenJets, b-tagging, grooming, njettiness)
process.load('BaconProd/Ntupler/myGenJets_cff')
process.load('BaconProd/Ntupler/myJetExtrasAK4CHS_cff')
process.load('BaconProd/Ntupler/myJetExtrasAK8CHS_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA8CHS_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA15CHS_cff')

process.load('BaconProd/Ntupler/myJetExtrasAK4Puppi_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA8Puppi_cff')
process.load('BaconProd/Ntupler/myJetExtrasCA15Puppi_cff')

from BaconProd.Ntupler.myGenJets_cff            import setMiniAODGenJets
from BaconProd.Ntupler.myJetExtrasAK4CHS_cff    import setMiniAODAK4CHS
from BaconProd.Ntupler.myJetExtrasAK8CHS_cff    import setMiniAODAK8CHS
from BaconProd.Ntupler.myJetExtrasCA8CHS_cff    import setMiniAODCA8CHS
from BaconProd.Ntupler.myJetExtrasCA15CHS_cff   import setMiniAODCA15CHS

from BaconProd.Ntupler.myJetExtrasAK4Puppi_cff  import setMiniAODAK4Puppi
from BaconProd.Ntupler.myJetExtrasCA8Puppi_cff  import setMiniAODCA8Puppi
from BaconProd.Ntupler.myJetExtrasCA15Puppi_cff import setMiniAODCA15Puppi

setMiniAODGenJets(process)
setMiniAODAK4CHS(process)
setMiniAODAK8CHS(process)
setMiniAODCA8CHS(process)
setMiniAODCA15CHS(process)

setMiniAODAK4Puppi (process)
setMiniAODCA8Puppi (process)
setMiniAODCA15Puppi(process)

# MVA MET
from BaconProd.Ntupler.myMVAMet_cff import setMiniAODMVAMet
process.load('BaconProd/Ntupler/myMVAMet_cff')     
setMiniAODMVAMet(process)

# PF MET corrections
process.load("BaconProd/Ntupler/myPFMETCorrections_cff")
process.pfJetMETcorr.jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3Corrector")
process.producePFMETCorrections = cms.Sequence(process.producePFMETCorrectionsMC)
if is_data_flag:
  process.pfJetMETcorr.jetCorrLabel = cms.InputTag("ak4PFL1FastL2L3ResidualCorrector")
  process.producePFMETCorrections = cms.Sequence(process.producePFMETCorrectionsData)
  process.AK4QGTaggerCHS.jec = cms.InputTag("ak4PFL1FastL2L3ResidualCorrector")

# PUPPI Woof Woof
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.puppi.candName       = cms.InputTag('packedPFCandidates')
process.puppi.vertexName     = cms.InputTag('offlineSlimmedPrimaryVertices')
process.pfCandNoLep = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("abs(pdgId) != 13 && abs(pdgId) != 11 && abs(pdgId) != 15"))
process.puppinolep = process.puppi.clone()
process.puppinolep.candName = 'pfCandNoLep'
process.load('RecoMET.METProducers.PFMET_cfi')
process.pfMetPuppi = process.pfMet.clone();
process.pfMetPuppi.src = cms.InputTag('puppinolep')
process.pfMetPuppi.calculateSignificance = False
process.pfJetMETcorrPuppi.jetCorrLabel = cms.InputTag("ak4PuppiL1FastL2L3Corrector")
process.producePFMETCorrectionsPuppi = cms.Sequence(process.producePFMETCorrectionsPuppiMC)
#MET 3.0
from BaconProd.Ntupler.myMET30_cff import setMet30
setMet30(process,3.0)

#--------------------------------------------------------------------------------
# input settings
#================================================================================
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )
process.source = cms.Source("PoolSource",
#fileNames = cms.untracked.vstring('/store/relval/CMSSW_7_4_6/RelValQCD_Pt_600_800_13/MINIAODSIM/MCRUN2_74_V9-v2/00000/9EFA69F6-E11A-E511-8002-00259059642A.root')#/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/70000/A2442E56-0409-E511-949E-0025905A60B0.root')
fileNames = cms.untracked.vstring('/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/00C4781D-6B08-E511-8A0A-0025905A6084.root','/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/00000/02DE3B74-6C08-E511-ABE3-0025905A60D0.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

#--------------------------------------------------------------------------------
# Reporting
#================================================================================
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)

#--------------------------------------------------------------------------------
# Bacon making settings
#================================================================================
process.ntupler = cms.EDAnalyzer('NtuplerMod',
  skipOnHLTFail     = cms.untracked.bool(do_hlt_filter),
  useAOD            = cms.untracked.bool(False),
  outputName        = cms.untracked.string('testbacon_miniaod.root'),
  TriggerFile       = cms.untracked.string(hlt_filename),
  edmPVName         = cms.untracked.string('offlineSlimmedPrimaryVertices'),
  edmGenRunInfoName = cms.untracked.string('generator'),
  
  Info = cms.untracked.PSet(
    isActive           = cms.untracked.bool(True),
    edmPFCandName      = cms.untracked.string('packedPFCandidates'),
    edmPileupInfoName  = cms.untracked.string('addPileupInfo'),
    edmBeamspotName    = cms.untracked.string('offlineBeamSpot'),
    edmMETName         = cms.untracked.string('slimmedMETs'),
    edmMVAMETName      = cms.untracked.string('pfMVAMEt'),
    edmPuppETName      = cms.untracked.string('pfMetPuppi'),
    edmPuppETCorrName  = cms.untracked.string('pfType1PuppiCorrectedMet'),
    edmPFMET30Name     = cms.untracked.string('pfMet30'),
    edmPFMETC30Name    = cms.untracked.string('pfType1CorrectedMet30'),
    edmMVAMET30Name    = cms.untracked.string('pfMVAMEt30'),
    edmPuppET30Name    = cms.untracked.string('pfMetPuppi30'),
    edmPuppET30CorrName = cms.untracked.string('pfType1PuppiCorrectedMet30'),
    edmRhoForIsoName   = cms.untracked.string('fixedGridRhoFastjetAll'),
    edmRhoForJetEnergy = cms.untracked.string('fixedGridRhoFastjetAll'),
    doFillMETFilters   = cms.untracked.bool(False),
    doFillMET          = cms.untracked.bool(True)
  ),
  
  GenInfo = cms.untracked.PSet(
    isActive            = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmGenEventInfoName = cms.untracked.string('generator'),
    edmGenParticlesName = cms.untracked.string('prunedGenParticles'),
    fillAllGen          = cms.untracked.bool(False),
    fillLHEWeights      = cms.untracked.bool(False)
  ),
  
  PV = cms.untracked.PSet(
    isActive      = cms.untracked.bool(True),   
    edmName       = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    minNTracksFit = cms.untracked.uint32(0),
    minNdof       = cms.untracked.double(4),
    maxAbsZ       = cms.untracked.double(24),
    maxRho        = cms.untracked.double(2)
  ),
  
  Electron = cms.untracked.PSet(
    isActive = cms.untracked.bool(True),
    minPt    = cms.untracked.double(7),
    edmName  = cms.untracked.string('slimmedElectrons'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppinolep'),
    usePuppi                  = cms.untracked.bool(True)
  ),
  
  Muon = cms.untracked.PSet(
    isActive = cms.untracked.bool(True),
    minPt    = cms.untracked.double(3),
    edmName  = cms.untracked.string('slimmedMuons'),
    #puppi
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppinolep'),
    usePuppi                  = cms.untracked.bool(True)    
  ),
  
  Photon = cms.untracked.PSet(
    isActive              = cms.untracked.bool(True),
    minPt                 = cms.untracked.double(10),
    edmName               = cms.untracked.string('slimmedPhotons'),
    edmChHadIsoMapTag     = cms.untracked.InputTag("photonIDValueMapProducer:phoChargedIsolation"),        # EGM recommendation not in AOD/MINIAOD
    edmNeuHadIsoMapTag    = cms.untracked.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),  # EGM recommendation not in AOD/MINIAOD
    edmGammaIsoMapTag     = cms.untracked.InputTag("photonIDValueMapProducer:phoPhotonIsolation")          # EGM recommendation not in AOD/MINIAOD
  ),
  
  Tau = cms.untracked.PSet(
    isActive = cms.untracked.bool(True),
    minPt    = cms.untracked.double(10),
    edmName  = cms.untracked.string('slimmedTaus'),
    edmPuppiName              = cms.untracked.string('puppi'),
    edmPuppiNoLepName         = cms.untracked.string('puppinolep'),
    usePuppi                  = cms.untracked.bool(True)
  ),
  
  AK4CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(False),
    minPt                = cms.untracked.double(20),
    coneSize             = cms.untracked.double(0.4),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmPVName            = cms.untracked.string('offlineSlimmedPrimaryVertices'),

    # names of various jet-related collections
    jetName              = cms.untracked.string('slimmedJets'),
    genJetName           = cms.untracked.string('slimmedGenJets'),
    csvBTagName          = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    qgLikelihood         = cms.untracked.string('QGTagger')
  ),

  AK4Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(20),
    coneSize             = cms.untracked.double(0.4),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    
    edmPVName = cms.untracked.string('offlineSlimmedPrimaryVertices'),
   # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK4puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK4puppi.txt',
				       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK4puppi.txt',
				       'BaconProd/Utils/data/MCRUN2_74_V9A_L2L3Residual_AK4puppi.txt')
                 if is_data_flag else 
                 cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK4puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK4puppi.txt',
				       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK4puppi.txt')
	       ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK4puppi.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK4puppi.txt')
                  ),
    edmRhoName = cms.untracked.string('fixedGridRhoFastjetAll'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    # names of various jet-related collections
    jetName            = cms.untracked.string('AK4PFJetsPuppi'),
    genJetName         = cms.untracked.string('AK4GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('AK4FlavorPuppi'),
    prunedJetName      = cms.untracked.string('AK4caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('AK4caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('AK4caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('AK4PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    jettiness          = cms.untracked.string('AK4NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('AK4QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('AK4QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('')
  ),

  AK8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(False),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmPVName            = cms.untracked.string('offlineSlimmedPrimaryVertices'),

    # names of various jet-related collections
    jetName              = cms.untracked.string('slimmedJetsAK8'),
    genJetName           = cms.untracked.string('AK8GenJetsCHS'),
    subJetName           = cms.untracked.string('SoftDrop'),
    csvBTagName          = cms.untracked.string('AK8PFCombinedInclusiveSecondaryVertexV2BJetTags'),
    qgLikelihood         = cms.untracked.string(''),
    prunedJetName        = cms.untracked.string('ak8PFJetsCHSPruned'),
    trimmedJetName       = cms.untracked.string('ak8PFJetsCHSTrimmed'),
    softdropJetName      = cms.untracked.string('ak8PFJetsCHSSoftDrop'),
    jettiness            = cms.untracked.string('NjettinessAK8'),
    topTaggerName        = cms.untracked.string('CMS')
  ),

  CA8CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(False),
    useAOD               = cms.untracked.bool(False),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(False),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmPVName            = cms.untracked.string('offlineSlimmedPrimaryVertices'),

    # names of various jet-related collections
    jetName              = cms.untracked.string('CA8PFJetsCHS'),
    genJetName           = cms.untracked.string('CA8GenJetsCHS'),
    subJetName           = cms.untracked.string('SoftDrop'),
    csvBTagName          = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTags'),
    qgLikelihood         = cms.untracked.string('QGTagger'),
    prunedJetName        = cms.untracked.string('CA8PFJetsCHSPruned'),
    trimmedJetName       = cms.untracked.string('CA8PFJetsCHSTrimmed'),
    softdropJetName      = cms.untracked.string('CA8PFJetsCHSSoftDrop'),
    jettiness            = cms.untracked.string('CA8NjettinessCHS'),
    topTaggerName        = cms.untracked.string('CMS')
  ),

 CA8Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(0.8),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    
    edmPVName = cms.untracked.string('offlineSlimmedPrimaryVertices'),
    
    # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2L3Residual_AK8puppi.txt')
                 if is_data_flag else
                 cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK8puppi.txt')
               ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK8puppi.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK8puppi.txt')
                  ),
    edmRhoName = cms.untracked.string('fixedGridRhoFastjetAll'),
    
    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),
    
    # names of various jet-related collections
    jetName            = cms.untracked.string('CA8PFJetsPuppi'),
    genJetName         = cms.untracked.string('CA8GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA8FlavorPuppi'),
    prunedJetName      = cms.untracked.string('CA8caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('CA8caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('CA8caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('CA8caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('CA8PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    jettiness          = cms.untracked.string('CA8NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('CA8QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('CA8QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('HEP')
  ),

  CA15CHS = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),
    edmPVName            = cms.untracked.string('offlineSlimmedPrimaryVertices'),

    # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK8PFchs.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK8PFchs.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK8PFchs.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2L3Residual_AK8PFchs.txt')
                 if is_data_flag else
                 cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK8PFchs.txt',   # (!) any JECs for CA15 CHS?
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK8PFchs.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK8PFchs.txt')
               ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK8PFchs.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK8PFchs.txt')
                  ),
    edmRhoName = cms.untracked.string('fixedGridRhoFastjetAll'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),

    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsCHS'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15FlavorCHS'),
    prunedJetName      = cms.untracked.string('CA15caPFJetsPrunedCHS'),
    trimmedJetName     = cms.untracked.string('CA15caPFJetsTrimmedCHS'),
    softdropJetName    = cms.untracked.string('CA15caPFJetsSoftDropCHS'),
    subJetName         = cms.untracked.string('CA15caPFJetsSoftDropCHS'),
    csvBTagName        = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsCHS'),
    csvBTagSubJetName  = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsSJCHS'),
    jettiness          = cms.untracked.string('CA15NjettinessCHS'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerCHS'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsCHS'),
    topTaggerName      = cms.untracked.string('HEP')
  ),
  CA15Puppi = cms.untracked.PSet(
    isActive             = cms.untracked.bool(True),
    useAOD               = cms.untracked.bool(True),
    minPt                = cms.untracked.double(180),
    coneSize             = cms.untracked.double(1.5),
    doComputeFullJetInfo = cms.untracked.bool(True),
    doGenJet             = ( cms.untracked.bool(False) if is_data_flag else cms.untracked.bool(True) ),

    edmPVName = cms.untracked.string('offlineSlimmedPrimaryVertices'),

    # ORDERED lists of jet energy correction input files
    jecFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK8Puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2L3Residual_AK8puppi.txt')
                 if is_data_flag else
                 cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_L1FastJet_AK8puppi.txt',   # (!) any JECs for CA15 Puppi?
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L2Relative_AK8puppi.txt',
                                       'BaconProd/Utils/data/MCRUN2_74_V9A_L3Absolute_AK8puppi.txt')
               ),
    jecUncFiles = ( cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK8puppi.txt')
                    if is_data_flag else
                    cms.untracked.vstring('BaconProd/Utils/data/MCRUN2_74_V9A_Uncertainty_AK8puppi.txt')
                  ),
    edmRhoName = cms.untracked.string('fixedGridRhoFastjetAll'),

    # ORDERD list of pileup jet ID input files
    jetPUIDFiles = cms.untracked.vstring('',
                                         'BaconProd/Utils/data/TMVAClassificationCategory_JetID_53X_chs_Dec2012.weights.xml'),

    # names of various jet-related collections
    jetName            = cms.untracked.string('CA15PFJetsPuppi'),
    genJetName         = cms.untracked.string('CA15GenJetsCHS'),
    jetFlavorName      = cms.untracked.string('CA15FlavorPuppi'),
    prunedJetName      = cms.untracked.string('CA15caPFJetsPrunedPuppi'),
    trimmedJetName     = cms.untracked.string('CA15caPFJetsTrimmedPuppi'),
    softdropJetName    = cms.untracked.string('CA15caPFJetsSoftDropPuppi'),
    subJetName         = cms.untracked.string('CA15caPFJetsSoftDropPuppi'),
    csvBTagName        = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsPuppi'),
    csvBTagSubJetName  = cms.untracked.string('CA15PFCombinedInclusiveSecondaryVertexV2BJetTagsSJPuppi'),
    jettiness          = cms.untracked.string('CA15NjettinessPuppi'),
    qgLikelihood       = cms.untracked.string('CA15QGTaggerPuppi'),
    qgLikelihoodSubjet = cms.untracked.string('CA15QGTaggerSubJetsPuppi'),
    topTaggerName      = cms.untracked.string('HEP')
  )
)

process.baconSequence = cms.Sequence(process.photonIDValueMapProducer *
                                     process.QGTagger                 *
                                     process.ak4PFJets                *
                                     process.packedPFCandidates30     *
                                     process.ak4PFJets30              *
                                     process.egmGsfElectronIDSequence * 
                                     process.egmPhotonIDSequence      * 
                                     process.slimmedMuonsTight        * 
                                     process.slimmedTausLoose         * 
                                     process.slimmedElectronsTight    * 
                                     process.pfMVAMEtSequenceNoLep    *
                                     process.pfMVAMEtSequenceNoLep30  *
                                     process.pfCandNoLep              *
                                     process.pfNoPileUpJME            *
                                     process.puppi                    *
                                     process.puppinolep               *
                                     process.puppi30                  *
                                     process.puppinolep30             *
                                     process.AK4PFJetsPuppi30         * 
                                     process.pfMetPuppi               *
                                     process.genjetsequence           *
                                     process.AK4genjetsequenceCHS     *
                                     #process.AK4jetsequenceCHS        *
                                     process.AK4jetsequencePuppi      *
                                     process.producePFMETCorrectionsPuppi*
                                     process.ak4PFL1FastL2L3CorrectorChain*
                                     process.allMET30                 *
                                     process.AK8jetsequenceCHS        *
                                     process.CA8jetsequenceCHS        *
                                     process.CA15jetsequenceCHS       *
                                     process.CA8jetsequencePuppi      *
                                     process.CA15jetsequencePuppi     *
                                     process.ntupler)

#--------------------------------------------------------------------------------
# apply trigger filter, if necessary
#================================================================================
if do_hlt_filter:
  process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')
  process.hltHighLevel.throw = cms.bool(False)
  process.hltHighLevel.HLTPaths = cms.vstring()
  hlt_file = open(cmssw_base + "/src/" + hlt_filename, "r")
  for line in hlt_file.readlines():
    line = line.strip()              # strip preceding and trailing whitespaces
    if (line[0:3] == 'HLT'):         # assumes typical lines begin with HLT path name (e.g. HLT_Mu15_v1)
      hlt_path = line.split()[0]
      process.hltHighLevel.HLTPaths.extend(cms.untracked.vstring(hlt_path))
  process.p = cms.Path(process.hltHighLevel*process.baconSequence)
else:
  process.p = cms.Path(process.baconSequence)

#--------------------------------------------------------------------------------
# simple checks to catch some mistakes...
#================================================================================
if is_data_flag:
  assert process.ntupler.GenInfo.isActive == cms.untracked.bool(False)
  assert process.ntupler.AK4CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA8CHS.doGenJet  == cms.untracked.bool(False)
  assert process.ntupler.CA15CHS.doGenJet == cms.untracked.bool(False)
