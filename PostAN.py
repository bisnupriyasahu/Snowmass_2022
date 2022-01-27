#!/usr/bin/env python

import sys

import ROOT
from array import array
import math
import re

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


def taucand1(tau,branchJet):
  i = 0
  tmpcand = []  
  #  tmpcand.clear()
  for tau in branchJet:
    i +=1
    #tmp_tau_2.SetPtEtaPhi(tau_pt[i],tau_eta[i],tau_phi[i]);
    wp = 2
    print "coming inside 1st tau"
    tautagOk = ( tau.TauTag & (1 << wp) )
    if (tau.PT >30 and abs(tau.Eta) < 3. and tautagOk):
      if (abs(tau.Charge) == 1):
        tmpcand = append(tau) 

  return tmpcand

def taucand2(tau,branchJet,index_tau1):
  i = 0
  tmpcand = []
  #tmpcand.clear()
  for tau in branchJet:
    i +=1
    if (index_tau1 >= 0 and tau == index_tau1): 
      continue
    #tmp_tau_2.SetPtEtaPhi(tau_pt[i],tau_eta[i],tau_phi[i]);                                                                                                                                               
    print "coming inside 2nd tau"
    wp = 2
    tautagOk = ( tau.TauTag & (1 << wp) )
    if (tau.PT >30 and abs(tau.Eta) < 3. and tautagOk):
      if (abs(tau.Charge) == 1):
        tmpcand = append(tau)
  return tmpcand



inputFile = sys.argv[1]
outputFile = sys.argv[2]
# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchJet      = treeReader.UseBranch("JetPUPPITight")
branchElectron = treeReader.UseBranch("ElectronMedium")
branchWeight   = treeReader.UseBranch("Weight")
branchEvent    = treeReader.UseBranch("Event")
branchPuppiJet = treeReader.UseBranch('JetPUPPI')
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

#histElectronPT = ROOT.TH1F("Electron_pt", "electron P_{T}", 100, 0.0, 1000.0)

# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  # If event contains at least 1 jet
  if branchJet.GetEntries() > 0:
    # Take first jet
    tau = branchJet.At(0)
    
    #opposite charged tau
    index_tau1 = index_tau1 = -1
    i = 0
    for tau in branchJet:
      i +=1
      index_tau1 = taucand1(tau,branchJet)
      print "index_tau 1  ", index_tau1.Charge
      index_tau2 = taucand2(tau,branchJet,index_tau1)  
      if( index_tau1 > -1 and index_tau2 > -1 ):
        if(tau.Charge(index_tau1) * tau.Charge(index_tau2) < 0):
          print "coming inside the taucharge"
'''
      if(tau_charge[index_tau1]*tau_charge[index_tau2] < 0  )
          
    ## 0 - Loose , 1 - Medium, 2 - Tight
    wp = 1

    BtagOk = ( tau.BTag & (1 << wp) )
    pt = tau.PT
    eta = abs(tau.Eta)

    # Plot jet transverse momentum
    if (BtagOk and pt > 30. and eta < 3.):
        tauPT.Fill(tau.PT)
'''       

#cnv.cd(2)
outputfile.cd()
tauPT.Write()
metPT.Write()
print tauPT.GetEntries()
print metPT.GetEntries()
outputfile.Close()
input("Press Enter to continue...")
