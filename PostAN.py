#!/Usr/bin/env python

import sys

import ROOT 
from array import array
import math
import re
from ROOT import *
import numpy as np
from mt2 import mt2

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

def nDaughters(gen):
  """Returns the number of daughters of a genparticle."""
  return gen.D2 - gen.D1

def finalDaughters(gen, brpart, daughters=None):
  if daughters is None:
    daughters = []
    
  for i in range(gen.D1, gen.D2+1):
    daughter = brpart[i]
    if nDaughters(daughter) == 0:
      daughters.append(daughter)
    else:
      finalDaughters(daughter, daughters)
    #print "daughter : ",daughter
  return daughters
      

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
nEvents = ROOT.TH1F("nEvents", "total events", 10, 0.0, 10)
tauPT_1 = ROOT.TH1F("tau1_pt", "tau P_{T}", 100, 30.0, 1000.0)
tauPT_2 = ROOT.TH1F("tau2_pt", "tau P_{T}", 100, 30.0, 1000.0)
metPT = ROOT.TH1F("MET", "MET", 50, 0.0, 500.0)
HT_Tot = ROOT.TH1F("HT", "Sum P_{T}", 100, 0.0, 1500.0)
ptratio_tau1 = ROOT.TH1F("ptratio_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
ptratio_tau2 = ROOT.TH1F("ptratio_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
genmatch_ptratio_tau1 = ROOT.TH1F("genmatch_ptratio_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
genmatch_ptratio_tau2 = ROOT.TH1F("genmatch_ptratio_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
notgenmatch_ptratio_tau1 = ROOT.TH1F("notgenmatch_ptratio_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
notgenmatch_ptratio_tau2 = ROOT.TH1F("notgenmatch_ptratio_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
MT = ROOT.TH1F("MT", "MT", 50, 0.0, 500.0)
DR_daughter = ROOT.TH1F("DeltaR","Delta R ",2000,0.0,3)


#histElectronPT = ROOT.TH1F("Electron_pt", "electron P_{T}", 100, 0.0, 1000.0)

# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)
  nEvents.Fill(1)
  nEvents.Fill(2,numberOfEntries)
  # If event contains at least 1 jet
  #if branchJet.GetEntries() > 0:
  # Take first jet
  tau1_idx = -1
  tau2_idx = -1
  tau1_tau2_HT = -1
  Tltau1_p4 = TLorentzVector()
  Tltau2_p4 = TLorentzVector()
  i = 0
  #print("coming in entries,", entry)
  for iTau1, tau1 in enumerate(branchJet) :
    i +=1
    #print("coming inside tau1")
    tautagOk1 = ( tau1.TauTag & (1 << 2) )
    if (not (tau1.PT >30 and abs(tau1.Eta) < 3. and tautagOk1 and abs(tau1.Charge) == 1)): continue
    for iTau2, tau2 in enumerate(branchJet):
      #print("coming inside tau2")
      if (iTau1 == iTau2) : continue
      tautagOk2 = ( tau2.TauTag & (1 << 2) )
      if (not (tau2.PT >30 and abs(tau2.Eta) < 3. and tautagOk2 and abs(tau2.Charge) == 1)): continue
      if (tau1.Charge*tau2.Charge < 0):
        HT = tau1.PT + tau2.PT
        if (HT > tau1_tau2_HT) :
          tau1_idx = iTau1
          tau2_idx = iTau2
          tau1_tau2_HT = HT
  
  if (not(tau1_idx >= 0 and tau2_idx >= 0)): continue
  tau1 = branchJet.At(tau1_idx)
  tau2 = branchJet.At(tau2_idx)
  if(tau1.PT < tau2.PT):
    tau = tau1
    tau1 = tau2
    tau2 = tau
  Tltau1_p4.SetPtEtaPhiM(tau1.PT, tau1.Eta, tau1.Phi, tau1.Mass)
  Tltau2_p4.SetPtEtaPhiM(tau2.PT, tau2.Eta, tau2.Phi, tau2.Mass)
  tau1tau2_m = (Tltau1_p4 + Tltau2_p4).M()

  btag_idx = -1    
  for ibjet, bjet in enumerate(branchJet) :
    if (ibjet == tau1_idx or ibjet == tau2_idx): continue
    
    btagok = (bjet.BTag & (1 << 1) )
    if (bjet.PT > 30 and abs(bjet.Eta) < 5. and btagok):
      btag_idx = ibjet
  #print( "btag index" , btag_idx)

  HT_Total = -1
  for ijet, jet in enumerate(branchJet) :
    HT_Total += jet.PT
  #print("Total HT", HT_Total)
    
  Met_PT = -1
  imet_idx = -1
  Met_Phi = 0
  for imet, met in enumerate(branchPuppiMissingET):
    Met_PT = met.MET
    Met_Phi = met.Phi
    imet_idx = imet
  #print(" imet_idx", imet_idx)
 

  
  if (not (tau1_idx >= 0 and tau2_idx >= 0 and btag_idx >= 0 and HT_Total > 100 and Met_PT > 50 and tau1tau2_m > 100)): continue
  print ("Invarient mass : ", tau1tau2_m)
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

  leadchtau1 = -1
  leadchtau2 = -1


  isLeptonic = False  
  tau1_leadCH = None
  all_consti1_p4 = TLorentzVector()
  for consti in tau1.Constituents:
    const_p4 =  TLorentzVector()
    const_p4.SetPtEtaPhiM(consti.PT, consti.Eta, consti.Phi, consti.Mass)
    all_consti1_p4 += const_p4
    ids = consti.PID
    if (abs(ids) in [11, 12, 13, 14]):
      isLeptonic = True
    if (isLeptonic == True or consti.Charge == 0) :
      continue
    const_p4 = TLorentzVector()
    const_p4.SetPtEtaPhiM(consti.PT, consti.Eta, consti.Phi, consti.Mass)
    dR = const_p4.DeltaR(Tltau1_p4)
    if (dR < 0.1):
      #print(dR)
      chpt = consti.PT
      if (tau1_leadCH is None or consti.PT > tau1_leadCH.PT):
        tau1_leadCH = consti
  if (tau1_leadCH is not None):
    leadchtau1 =  tau1_leadCH.PT/all_consti1_p4.Pt()
    ptratio_tau1.Fill(leadchtau1)

  isLeptonic2 = False  
  tau2_leadCH = None
  all_consti2_p4= TLorentzVector()
  for consti2 in tau2.Constituents:
    const2_p4 = TLorentzVector()
    const2_p4.SetPtEtaPhiM(consti2.PT, consti2.Eta, consti2.Phi, consti2.Mass)
    all_consti2_p4 += const2_p4
    ids2 = consti2.PID
    if (abs(ids2) in [11, 12, 13, 14]):
      isLeptonic2 = True
    if (isLeptonic2 == True or consti2.Charge == 0) :
      continue
    const2_p4 = TLorentzVector()
    const2_p4.SetPtEtaPhiM(consti2.PT, consti2.Eta, consti2.Phi, consti2.Mass)
    dR2 = const2_p4.DeltaR(Tltau2_p4)
    
    if (dR2 < 0.1):
      #print("coming inside 3 dr2 :",dR2)  
      chpt = consti.PT
      if (tau2_leadCH is None or consti2.PT > tau2_leadCH.PT):
        tau2_leadCH = consti2
  if (tau2_leadCH is not None):
    leadchtau2 =  tau2_leadCH.PT/all_consti2_p4.Pt()
    #print("leadchtau2 is ",leadchtau2)
    ptratio_tau2.Fill(leadchtau2)
  

 
  gen_1 = None
  gen_2 = None
  dr_dau = -1
  #  gen_elep4 = gen_mup4 = gen_taup4 = TLorentzVector()
  gen1_p4 = gen2_p4 = TLorentzVector()
  nele = 0
  #print("Branch part : ", branchParticle.GetEntries(), entry)
  #print( "Branch of jets :", branchJet.GetEntries(), entry)
  for igen,gen in enumerate(branchParticle):
    if(abs(gen.PID) == 15):
      gen_p4 = TLorentzVector()
      gen_p4.SetPtEtaPhiM(gen.PT, gen.Eta, gen.Phi, gen.Mass)
      dr_1 = gen_p4.DeltaR(Tltau1_p4)
      dr_2 = gen_p4.DeltaR(Tltau2_p4)
      if (dr_1 < 0.3):
        gen_1 = gen
      
      elif (dr_2 < 0.3):
        gen_2 = gen
          
  if (gen_1 is not None):
    #print (leadchtau1)
    gen_1pt =  gen_1.PT/tau1.PT
    gen1_p4.SetPtEtaPhiM(gen_1.PT, gen_1.Eta, gen_1.Phi, gen_1.Mass)     

    genmatch_ptratio_tau1.Fill(leadchtau1)
  else: 
    print (leadchtau1)
    notgenmatch_ptratio_tau1.Fill(leadchtau1)
  if (gen_2 is not None):
    gen_2pt =  gen_2.PT/tau2.PT
    gen2_p4.SetPtEtaPhiM(gen_2.PT, gen_2.Eta, gen_2.Phi, gen_2.Mass)     

    genmatch_ptratio_tau2.Fill(leadchtau2)
  else:
    print (leadchtau2)
    notgenmatch_ptratio_tau2.Fill(leadchtau2)




  
  lep_p4 = TLorentzVector()
  lep1_dr = lep2_dr = -1
  min_dr = 999.9 
  for jgen,gen in enumerate(branchParticle):

    if(abs(gen.PID) in [11,13]):
      lep_p4.SetPtEtaPhiM(gen.PT, gen.Eta, gen.Phi, gen.Mass)
      if (gen_1 is not None): 
        lep1_dr = lep_p4.DeltaR(gen1_p4)       
      elif(gen_2 is not None):   
        lep2_dr = lep_p4.DeltaR(gen2_p4)
      if (lep1_dr < min_dr):
        min_dr = lep1_dr
      elif (lep2_dr < min_dr):
        min_dr = lep2_dr
  DR_daughter.Fill(min_dr)
  print("min dr is : ",min_dr)


  tau1_px = Tltau1_p4.Px()
  tau1_py = Tltau1_p4.Py()
  tau1_m = tau_1.Mass
  #print("tau1_px and tau1_py", tau1_px, tau1_py)

  tau2_px = Tltau2_p4.Px()
  tau2_py = Tltau2_p4.Py()
  tau2_m = tau_2.Mass
  #print("tau2_px and tau2_py", tau2_px, tau2_py)

  MET_px = Met_PT*(math.cos(Met_Phi))
  MET_py = Met_PT*(math.sin(Met_Phi))
  #print("Met_px and MET_py", MET_px, MET_py)


  MT_ = mt2(
    tau1_m,tau1_px,tau1_py,
    tau2_m,tau2_px,tau2_py,
    MET_px,MET_py,
    0,0)
  print("MT is : ", MT_)
  MT.Fill(MT_)


#cnv.cd(2)
outputfile.cd()
nEvents.Write()
tauPT_1.Write()
tauPT_2.Write()
metPT.Write()
HT_Tot.Write()
ptratio_tau1.Write()
ptratio_tau2.Write()
genmatch_ptratio_tau1.Write()
genmatch_ptratio_tau2.Write()
notgenmatch_ptratio_tau1.Write()
notgenmatch_ptratio_tau2.Write()
MT.Write()
DR_daughter.Write()

print( nEvents.GetEntries())
print( tauPT_1.GetEntries())
print( tauPT_2.GetEntries())
print( metPT.GetEntries())
print( HT_Tot.GetEntries())
print( ptratio_tau1.GetEntries())
print( ptratio_tau2.GetEntries())
print( genmatch_ptratio_tau1.GetEntries())
print( genmatch_ptratio_tau2.GetEntries())
print( notgenmatch_ptratio_tau1.GetEntries())
print( notgenmatch_ptratio_tau2.GetEntries())
print( MT.GetEntries())
print( DR_daughter.GetEntries())

outputfile.Close()
input("Press Enter to continue...")

