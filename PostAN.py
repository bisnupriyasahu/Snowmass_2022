#!/usr/bin/env python

import sys

import ROOT
from array import array
import math
import re

import numpy as np


try:
  input = raw_input
except:
  pass

if len(sys.argv) < 2:
  print(" Usage: Example1.py input_file")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass

inputFile = sys.argv[1]
outputFile = sys.argv[2]
# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
#branchJet      = treeReader.UseBranch("JetPUPPITight")
branchElectron = treeReader.UseBranch("ElectronMedium")
branchWeight   = treeReader.UseBranch("Weight")
branchEvent    = treeReader.UseBranch("Event")
branchJet = treeReader.UseBranch('JetPUPPI')
branchPuppiMissingET  = treeReader.UseBranch('PuppiMissingET')
branchPuppiCandidate  = treeReader.UseBranch('ParticleFlowCandidate')
branchRho             = treeReader.UseBranch('Rho')
branchParticle        = treeReader.UseBranch('Particle') 
branchGenJet          = treeReader.UseBranch('GenJet')
branchPuppiJetLoose   = treeReader.UseBranch('JetPUPPILoose')
#branchPuppiJetTight   = treeReader.UseBranch('JetPUPPITight')
branchFatJet          = treeReader.UseBranch('JetPUPPIAK8')


# Book histograms
outputfile = ROOT.TFile(outputFile, 'RECREATE')
tauPT = ROOT.TH1F("tau_pt", "tau P_{T}", 100, 0.0, 1000.0)
metPT = ROOT.TH1F("met_pt", "met P_{T}", 100, 0.0, 1000.0)
HT_Tot = ROOT.TH1F("HT", "P_{T}", 100, 0.0, 1000.0)
#histElectronPT = ROOT.TH1F("Electron_pt", "electron P_{T}", 100, 0.0, 1000.0)

# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  # If event contains at least 1 jet
  #if branchJet.GetEntries() > 0:
  # Take first jet
  tau1_idx = -1
  tau2_idx = -1
  tau1_tau2_HT = -1

  i = 0
  for iTau1, tau1 in enumerate(branchJet) :
    i +=1
    tautagOk1 = ( tau1.TauTag & (1 << 2) )
    if (tau1.PT >30 and abs(tau1.Eta) < 3. and tautagOk1):
      if (abs(tau1.Charge) != 1): continue
      #else: print "tau1.Charge", tau1.Charge
      for iTau2, tau2 in enumerate(branchJet) :
        if (iTau1 == iTau2) : continue
        tautagOk2 = ( tau2.TauTag & (1 << 2) )
        if (tau2.PT >30 and abs(tau2.Eta) < 3. and tautagOk2):
          if (abs(tau2.Charge) != 1):
            continue
          #else: print "tau2.Charge", tau2.Charge
      if (tau1.Charge*tau2.Charge < 0):
        HT = tau1.PT + tau2.PT
        print "HT in loop ", HT
        if (HT > tau1_tau2_HT) :
            tau1_idx = iTau1
            tau2_idx = iTau2
            tau1_tau2_HT = HT
  btag_idx = -1    
  for ibjet, bjet in enumerate(branchJet) :
    if (ibjet == tau1_idx or ibjet == tau2_idx): continue
    
    btagok = (bjet.BTag & (1 << 1) )
    if (bjet.PT > 30 and abs(bjet.Eta) < 5. and btagok):
      btag_idx = ibjet
  print "btag index" , btag_idx

  HT_Total = -1
  for ijet, jet in enumerate(branchJet) :
    HT_Total += jet.PT
  print "Total HT", HT_Total
    
  Met_PT = -1
  imet_idx = -1
  for imet, met in enumerate(branchPuppiMissingET):
    if (met.MET > 50):
      Met_PT = met.MET
      imet_idx = imet

  if (not (tau1_idx > 0 and tau2_idx and btag_idx > 0 and HT_Total > 100 and Met_PT > 50)): continue
  for j, tau in enumerate(branchJet) :
    if (j == tau1_idx or j == tau2_idx or j == btag_idx):
      tauPT.Fill(tau.PT)
      metPT.Fill(Met_PT)
      HT_Tot.Fill(HT_Total)
  





#cnv.cd(2)
outputfile.cd()
tauPT.Write()
metPT.Write()
HT_Tot.Write()
print tauPT.GetEntries()
print metPT.GetEntries()
print HT.GetEntries()

outputfile.Close()
input("Press Enter to continue...")
