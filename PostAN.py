#!/Usr/bin/env python

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
tauPT_1 = ROOT.TH1F("tau1_pt", "tau P_{T}", 100, 0.0, 1000.0)
tauPT_2 = ROOT.TH1F("tau2_pt", "tau P_{T}", 100, 0.0, 1000.0)
metPT = ROOT.TH1F("met_pt", "met P_{T}", 100, 0.0, 1000.0)
HT_Tot = ROOT.TH1F("HT", "P_{T}", 100, 0.0, 1000.0)
ptratio_tau1 = ROOT.TH1F("ptratio_tau1", "P_{T}", 50, 0.0, 2.0)
ptratio_tau2 = ROOT.TH1F("ptratio_tau2", "P_{T}", 50, 0.0, 2.0)
gen_ptratio_tau1 = ROOT.TH1F("gen_ptratio_tau1", "P_{T}", 50, 0.0, 2.0)
gen_ptratio_tau2 = ROOT.TH1F("gen_ptratio_tau2", "P_{T}", 50, 0.0, 2.0)




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
        #print "HT in loop ", HT
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
  #print "btag index" , btag_idx

  HT_Total = -1
  for ijet, jet in enumerate(branchJet) :
    HT_Total += jet.PT
  #print "Total HT", HT_Total
    
  Met_PT = -1
  imet_idx = -1
  for imet, met in enumerate(branchPuppiMissingET):
    Met_PT = met.MET
    imet_idx = imet
  #print imet_idx
  

  if (not (tau1_idx > 0 and tau2_idx and btag_idx > 0 and HT_Total > 100 and Met_PT > 50)): continue
  tau_1 = branchJet.At(tau1_idx)
  tau_2 = branchJet.At(tau2_idx)
  met_pt = branchPuppiMissingET.At(imet_idx)
  tau1pt = tau_1.PT
  tau2pt = tau_2.PT
  metpt = met_pt.MET

  tauPT_1.Fill(tau1pt)
  tauPT_2.Fill(tau2pt)
  metPT.Fill(metpt)
  HT_Tot.Fill(HT_Total)


  isLeptonic = False  
  #CH_PT = -1
  tau1_leadCH = None
  for consti in tau1.Constituents:
    ids = consti.PID
    #CH_PT = consti.PT
    if (abs(ids) in [11, 12, 13, 14]):
      isLeptonic = True
   
    if (isLeptonic == True or consti.Charge == 0) :
      continue
    #print "ids of tau1", consti.PID
    dR = math.sqrt(  (consti.Eta-(tau1.Eta)) * ((consti.Eta)-(tau1.Eta)) + ((consti.Phi)-(tau1.Phi)) * ((consti.Phi)-(tau1.Phi)) )
    if (dR < 0.1):
      print dR
      chpt = consti.PT
      #print "all ch pt is ", chpt
      if (tau1_leadCH is None or consti.PT > tau1_leadCH.PT):
        tau1_leadCH = consti
  if (tau1_leadCH is not None):
    leadchtau1 =  tau1_leadCH.PT/tau1.PT
    ptratio_tau1.Fill(leadchtau1)

  isLeptonic2 = False  
  #CH_PT = -1
  tau2_leadCH = None
  for consti2 in tau1.Constituents:
    ids2 = consti2.PID
    #CH_PT = consti.PT
    if (abs(ids2) in [11, 12, 13, 14]):
      isLeptonic2 = True
   
    if (isLeptonic2 == True or consti2.Charge == 0) :
      continue
    #print "ids of tau1", consti.PID
    dR2 = math.sqrt(  (consti.Eta-(tau1.Eta)) * ((consti.Eta)-(tau1.Eta)) + ((consti.Phi)-(tau1.Phi)) * ((consti.Phi)-(tau1.Phi)) )
    if (dR2 < 0.1):
      #print dR
      chpt = consti.PT
      #print "all ch pt is ", chpt
      if (tau2_leadCH is None or consti2.PT > tau2_leadCH.PT):
        tau2_leadCH = consti2
  if (tau2_leadCH is not None):
    leadchtau2 =  tau2_leadCH.PT/tau2.PT
    ptratio_tau2.Fill(leadchtau2)
  

 
  gen_1 = None
  gen_2 = None
  for igen,gen in enumerate(branchParticle):
    if (abs(gen.PID) == 15):
      print "gen pid ",gen.PID
      dr_1 =  math.sqrt(  (gen.Eta-(tau1.Eta)) * ((gen.Eta)-(tau1.Eta)) + ((gen.Phi)-(tau1.Phi)) * ((gen.Phi)-(tau1.Phi)) )
      dr_2 =  math.sqrt(  (gen.Eta-(tau2.Eta)) * ((gen.Eta)-(tau2.Eta)) + ((gen.Phi)-(tau2.Phi)) * ((gen.Phi)-(tau2.Phi)) )
      print "dr 1 is ", dr_1
      print "dr 2 is ",dr_2

      if (dr_1 < 0.1):
        gen_1 = gen
      if (dr_2 < 0.1):
        gen_2 = gen
        print "gen2 pt ", gen.PT
  if (gen_1 is not None):
    gen_1pt =  gen_1.PT/tau1.PT
    gen_ptratio_tau1.Fill(gen_1pt)
  if (gen_2 is not None):
    gen_2pt =  gen_2.PT/tau2.PT
    gen_ptratio_tau2.Fill(gen_2pt)





#cnv.cd(2)
outputfile.cd()
tauPT_1.Write()
tauPT_2.Write()
metPT.Write()
HT_Tot.Write()
ptratio_tau1.Write()
ptratio_tau2.Write()
gen_ptratio_tau1.Write()
gen_ptratio_tau2.Write()

print tauPT_1.GetEntries()
print tauPT_2.GetEntries()
print metPT.GetEntries()
print HT_Tot.GetEntries()
print ptratio_tau1.GetEntries()
print ptratio_tau2.GetEntries()
print gen_ptratio_tau1.GetEntries()
print gen_ptratio_tau2.GetEntries()

outputfile.Close()
input("Press Enter to continue...")
