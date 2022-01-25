//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 17 07:49:45 2022 by ROOT version 6.24/06
// from TTree mytree/TestTree
// found on file: TTJets_DiLept_TuneCUETP8M1_14TeV-madgraphMLM-pythia8_ntuple_93_0.root
//////////////////////////////////////////////////////////

#ifndef TauClass_h
#define TauClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.                                                                                                                                                 
#include "vector"
#include <vector>
#include <iostream>
#include <fstream>
#include <TH2.h>
#include <TH1.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TChain.h>

using namespace std;


class TauClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TFile *fileName;
   TTree *tree;

   //Histograms declaration
   TH1F *tau_pt_;
   TH1F *tau_eta_;
   TH1F *met_pt_;
   TH1F *met_phi_;
   TH1F *tauVis_pt_;
   TH1F *tauVis_M_;   
   TH1F *tauVis_eta_;
   TH1F *pi_pt_;
   TH1F *ptratio_;
   TH1F *pipt_tauvispt_;
   TH1F *tauvispt_taupt_;
   TH1F *HT;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evt_size;
   Int_t           vtx_size;

   Int_t           PFcandidate_size;
   Float_t         PFcandidate_pid[416];   //[PFcandidate_size]                                                                                                                                          
   Float_t         PFcandidate_pt[416];   //[PFcandidate_size]                                                                                                                                           
   Float_t         PFcandidate_eta[416];   //[PFcandidate_size]                                                                                                                                          
   Float_t         PFcandidate_phi[416];   //[PFcandidate_size]                                                                                                                                          
   Float_t         PFcandidate_mass[416];   //[PFcandidate_size]                                                                                                                                         
   Float_t         PFcandidate_charge[416];   //[PFcandidate_size]   
   Int_t           trueInteractions;
   Int_t           npuVertices;
   Float_t         vtx_pt2[258];   //[vtx_size]
   Float_t         vtx_x[258];   //[vtx_size]
   Float_t         vtx_y[258];   //[vtx_size]
   Float_t         vtx_z[258];   //[vtx_size]
   Float_t         genweight;
   Int_t           lheweight_size;
   Int_t           lheweight_val[446];   //[lheweight_size]
   Int_t           genpart_size;
   Int_t           genpart_pid[612];   //[genpart_size]
   Int_t           genpart_status[612];   //[genpart_size]
   Float_t         genpart_pt[612];   //[genpart_size]
   Float_t         genpart_eta[612];   //[genpart_size]
   Float_t         genpart_phi[612];   //[genpart_size]
   Float_t         genpart_mass[612];   //[genpart_size]
   Int_t           genpart_m1[612];   //[genpart_size]
   Int_t           genpart_m2[612];   //[genpart_size]
   Int_t           genpart_d1[612];   //[genpart_size]
   Int_t           genpart_d2[612];   //[genpart_size]
   Int_t           genjet_size;
   Float_t         genjet_pt[17];   //[genjet_size]
   Float_t         genjet_eta[17];   //[genjet_size]
   Float_t         genjet_phi[17];   //[genjet_size]
   Float_t         genjet_mass[17];   //[genjet_size]
   Int_t           genmet_size;
   Float_t         genmet_pt[1];   //[genmet_size]
   Float_t         genmet_phi[1];   //[genmet_size]
   Int_t           gamma_size;
   Float_t         gamma_pt[41];   //[gamma_size]
   Float_t         gamma_eta[41];   //[gamma_size]
   Float_t         gamma_phi[41];   //[gamma_size]
   Float_t         gamma_mass[41];   //[gamma_size]
   Float_t         gamma_idvar[41];   //[gamma_size]
   Float_t         gamma_reliso[41];   //[gamma_size]
   UInt_t          gamma_idpass[41];   //[gamma_size]
   UInt_t          gamma_isopass[41];   //[gamma_size]
   Int_t           elec_size;
   Float_t         elec_pt[31];   //[elec_size]
   Float_t         elec_eta[31];   //[elec_size]
   Float_t         elec_phi[31];   //[elec_size]
   Float_t         elec_mass[31];   //[elec_size]
   Int_t           elec_charge[31];   //[elec_size]
   Float_t         elec_idvar[31];   //[elec_size]
   Float_t         elec_reliso[31];   //[elec_size]
   UInt_t          elec_idpass[31];   //[elec_size]
   UInt_t          elec_isopass[31];   //[elec_size]
   Int_t           muon_size;
   Float_t         muon_pt[16];   //[muon_size]
   Float_t         muon_eta[16];   //[muon_size]
   Float_t         muon_phi[16];   //[muon_size]
   Float_t         muon_mass[16];   //[muon_size]
   Int_t           muon_charge[16];   //[muon_size]
   Float_t         muon_idvar[16];   //[muon_size]
   Float_t         muon_reliso[16];   //[muon_size]
   UInt_t          muon_idpass[16];   //[muon_size]
   UInt_t          muon_isopass[16];   //[muon_size]
   Int_t           tau_size;
   Float_t         tau_pt[14];   //[tau_size]
   Float_t         tau_eta[14];   //[tau_size]
   Float_t         tau_phi[14];   //[tau_size]
   Float_t         tau_mass[14];   //[tau_size]
   Int_t           tau_charge[14];   //[tau_size]
   Float_t         tau_decaymode[14];   //[tau_size]
   Float_t         tau_neutraliso[14];   //[tau_size]
   Float_t         tau_chargediso[14];   //[tau_size]
   Float_t         tau_combinediso[14];   //[tau_size]
   UInt_t          tau_isopass[14];   //[tau_size]
   Int_t           jetpuppi_size;
   Float_t         jetpuppi_pt[14];   //[jetpuppi_size]
   Float_t         jetpuppi_eta[14];   //[jetpuppi_size]
   Float_t         jetpuppi_phi[14];   //[jetpuppi_size]
   Float_t         jetpuppi_mass[14];   //[jetpuppi_size]
   UInt_t          jetpuppi_idpass[14];   //[jetpuppi_size]
   Float_t         jetpuppi_DeepJET[14];   //[jetpuppi_size]
   Int_t           jetpuppi_btag[14];   //[jetpuppi_size]
   Int_t           jetchs_size;
   Float_t         jetchs_pt[1];   //[jetchs_size]
   Float_t         jetchs_eta[1];   //[jetchs_size]
   Float_t         jetchs_phi[1];   //[jetchs_size]
   Float_t         jetchs_mass[1];   //[jetchs_size]
   UInt_t          jetchs_idpass[1];   //[jetchs_size]
   Float_t         jetchs_DeepJET[1];   //[jetchs_size]
   Int_t           jetchs_btag[1];   //[jetchs_size]
   Int_t           fatjet_size;
   //   Float_t         fatjet_pt[4];   //[fatjet_size]
   Float_t         fatjet_eta[4];   //[fatjet_size]
   Float_t         fatjet_phi[4];   //[fatjet_size]
   Float_t         fatjet_mass[4];   //[fatjet_size]
   Float_t         fatjet_pt[4];   //[fatjet_size]
   Float_t         fatjet_tau1[4];   //[fatjet_size]
   Float_t         fatjet_tau2[4];   //[fatjet_size]
   Float_t         fatjet_tau3[4];   //[fatjet_size]
   Float_t         fatjet_tau4[4];   //[fatjet_size]
   Float_t         fatjet_msoftdrop[4];   //[fatjet_size]
   Int_t           metpuppi_size;
   Float_t         metpuppi_pt[1];   //[metpuppi_size]
   Float_t         metpuppi_phi[1];   //[metpuppi_size]
   Int_t           metpf_size;
   Float_t         metpf_pt[1];   //[metpf_size]
   Float_t         metpf_phi[1];   //[metpf_size]

   // List of branches
   TBranch        *b_evt_size;   //!
   TBranch        *b_vtx_size;   //!
   TBranch        *b_PFcandidate_size;   //!                                                                                                                                                             
   TBranch        *b_PFcandidate_pid;   //!                                                                                                                                                              
   TBranch        *b_PFcandidate_pt;   //!                                                                                                                                                               
   TBranch        *b_PFcandidate_eta;   //!                                                                                                                                                              
   TBranch        *b_PFcandidate_phi;   //!                                                                                                                                                              
   TBranch        *b_PFcandidate_mass;   //!                                                                                                                                                             
   TBranch        *b_PFcandidate_charge;   //!         
   TBranch        *b_trueInteractions;   //!
   TBranch        *b_npuVertices;   //!
   TBranch        *b_vtx_pt2;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_genweight;   //!
   TBranch        *b_lheweight_size;   //!
   TBranch        *b_lheweight_val;   //!
   TBranch        *b_genpart_size;   //!
   TBranch        *b_genpart_pid;   //!
   TBranch        *b_genpart_status;   //!
   TBranch        *b_genpart_pt;   //!
   TBranch        *b_genpart_eta;   //!
   TBranch        *b_genpart_phi;   //!
   TBranch        *b_genpart_mass;   //!
   TBranch        *b_genpart_m1;   //!
   TBranch        *b_genpart_m2;   //!
   TBranch        *b_genpart_d1;   //!
   TBranch        *b_genpart_d2;   //!
   TBranch        *b_genjet_size;   //!
   TBranch        *b_genjet_pt;   //!
   TBranch        *b_genjet_eta;   //!
   TBranch        *b_genjet_phi;   //!
   TBranch        *b_genjet_mass;   //!
   TBranch        *b_genmet_size;   //!
   TBranch        *b_genmet_pt;   //!
   TBranch        *b_genmet_phi;   //!
   TBranch        *b_gamma_size;   //!
   TBranch        *b_gamma_pt;   //!
   TBranch        *b_gamma_eta;   //!
   TBranch        *b_gamma_phi;   //!
   TBranch        *b_gamma_mass;   //!
   TBranch        *b_gamma_idvar;   //!
   TBranch        *b_gamma_reliso;   //!
   TBranch        *b_gamma_idpass;   //!
   TBranch        *b_gamma_isopass;   //!
   TBranch        *b_elec_size;   //!
   TBranch        *b_elec_pt;   //!
   TBranch        *b_elec_eta;   //!
   TBranch        *b_elec_phi;   //!
   TBranch        *b_elec_mass;   //!
   TBranch        *b_elec_charge;   //!
   TBranch        *b_elec_idvar;   //!
   TBranch        *b_elec_reliso;   //!
   TBranch        *b_elec_idpass;   //!
   TBranch        *b_elec_isopass;   //!
   TBranch        *b_muon_size;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_mass;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_idvar;   //!
   TBranch        *b_muon_reliso;   //!
   TBranch        *b_muon_idpass;   //!
   TBranch        *b_muon_isopass;   //!
   TBranch        *b_tau_size;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_tau_eta;   //!
   TBranch        *b_tau_phi;   //!
   TBranch        *b_tau_mass;   //!
   TBranch        *b_tau_charge;   //!
   TBranch        *b_tau_decaymode;   //!
   TBranch        *b_tau_neutraliso;   //!
   TBranch        *b_tau_chargediso;   //!
   TBranch        *b_tau_combinediso;   //!
   TBranch        *b_tau_isopass;   //!
   TBranch        *b_jetpuppi_size;   //!
   TBranch        *b_jetpuppi_pt;   //!
   TBranch        *b_jetpuppi_eta;   //!
   TBranch        *b_jetpuppi_phi;   //!
   TBranch        *b_jetpuppi_mass;   //!
   TBranch        *b_jetpuppi_idpass;   //!
   TBranch        *b_jetpuppi_DeepJET;   //!
   TBranch        *b_jetpuppi_btag;   //!
   TBranch        *b_jetchs_size;   //!
   TBranch        *b_jetchs_pt;   //!
   TBranch        *b_jetchs_eta;   //!
   TBranch        *b_jetchs_phi;   //!
   TBranch        *b_jetchs_mass;   //!
   TBranch        *b_jetchs_idpass;   //!
   TBranch        *b_jetchs_DeepJET;   //!
   TBranch        *b_jetchs_btag;   //!
   TBranch        *b_fatjet_size;   //!
   //TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_eta;   //!
   TBranch        *b_fatjet_phi;   //!
   TBranch        *b_fatjet_mass;   //!
   TBranch        *b_fatjet_pt;   //!
   TBranch        *b_fatjet_tau1;   //!
   TBranch        *b_fatjet_tau2;   //!
   TBranch        *b_fatjet_tau3;   //!
   TBranch        *b_fatjet_tau4;   //!
   TBranch        *b_fatjet_msoftdrop;   //!
   TBranch        *b_metpuppi_size;   //!
   TBranch        *b_metpuppi_pt;   //!
   TBranch        *b_metpuppi_phi;   //!
   TBranch        *b_metpf_size;   //!
   TBranch        *b_metpf_pt;   //!
   TBranch        *b_metpf_phi;   //!

   TauClass(const char* file1, const char* file2);
   virtual ~TauClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Long64_t maxEvents, int reportEvery);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void Histos(const char* file2);
   virtual int TauCand1();
   virtual int TauCand2(int tau1Index);
   virtual vector<int> BTagjet();
   


};

#endif

#ifdef TauClass_cxx
TauClass::TauClass(const char* file1, const char* file2)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   

  TChain *chain = new TChain("myana/mytree");
  TString path = file1;


  TSystemDirectory sourceDir("hi",path);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter next(fileList);
  TSystemFile* filename;
  int fileNumber = 0;
  int maxFiles = -1;
  std::cout<<"path:"<<path<<std::endl;

  chain->Add(file1);
  std::cout<<"events inside .h"<<chain->GetEntry()<<std::endl;

  Init(chain);
  Histos(file2);

  }

TauClass::~TauClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   fileName->Close();
}

Int_t TauClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TauClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TauClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt_size", &evt_size, &b_evt_size);
   fChain->SetBranchAddress("vtx_size", &vtx_size, &b_vtx_size);
   fChain->SetBranchAddress("PFcandidate_size", &PFcandidate_size, &b_PFcandidate_size);
   fChain->SetBranchAddress("PFcandidate_pid", PFcandidate_pid, &b_PFcandidate_pid);
   fChain->SetBranchAddress("PFcandidate_pt", PFcandidate_pt, &b_PFcandidate_pt);
   fChain->SetBranchAddress("PFcandidate_eta", PFcandidate_eta, &b_PFcandidate_eta);
   fChain->SetBranchAddress("PFcandidate_phi", PFcandidate_phi, &b_PFcandidate_phi);
   fChain->SetBranchAddress("PFcandidate_mass", PFcandidate_mass, &b_PFcandidate_mass);
   fChain->SetBranchAddress("PFcandidate_charge", PFcandidate_charge, &b_PFcandidate_charge);
   fChain->SetBranchAddress("trueInteractions", &trueInteractions, &b_trueInteractions);
   fChain->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
   fChain->SetBranchAddress("vtx_pt2", vtx_pt2, &b_vtx_pt2);
   fChain->SetBranchAddress("vtx_x", vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("genweight", &genweight, &b_genweight);
   fChain->SetBranchAddress("lheweight_size", &lheweight_size, &b_lheweight_size);
   fChain->SetBranchAddress("lheweight_val", lheweight_val, &b_lheweight_val);
   fChain->SetBranchAddress("genpart_size", &genpart_size, &b_genpart_size);
   fChain->SetBranchAddress("genpart_pid", genpart_pid, &b_genpart_pid);
   fChain->SetBranchAddress("genpart_status", genpart_status, &b_genpart_status);
   fChain->SetBranchAddress("genpart_pt", genpart_pt, &b_genpart_pt);
   fChain->SetBranchAddress("genpart_eta", genpart_eta, &b_genpart_eta);
   fChain->SetBranchAddress("genpart_phi", genpart_phi, &b_genpart_phi);
   fChain->SetBranchAddress("genpart_mass", genpart_mass, &b_genpart_mass);
   fChain->SetBranchAddress("genpart_m1", genpart_m1, &b_genpart_m1);
   fChain->SetBranchAddress("genpart_m2", genpart_m2, &b_genpart_m2);
   fChain->SetBranchAddress("genpart_d1", genpart_d1, &b_genpart_d1);
   fChain->SetBranchAddress("genpart_d2", genpart_d2, &b_genpart_d2);
   fChain->SetBranchAddress("genjet_size", &genjet_size, &b_genjet_size);
   fChain->SetBranchAddress("genjet_pt", genjet_pt, &b_genjet_pt);
   fChain->SetBranchAddress("genjet_eta", genjet_eta, &b_genjet_eta);
   fChain->SetBranchAddress("genjet_phi", genjet_phi, &b_genjet_phi);
   fChain->SetBranchAddress("genjet_mass", genjet_mass, &b_genjet_mass);
   fChain->SetBranchAddress("genmet_size", &genmet_size, &b_genmet_size);
   fChain->SetBranchAddress("genmet_pt", genmet_pt, &b_genmet_pt);
   fChain->SetBranchAddress("genmet_phi", genmet_phi, &b_genmet_phi);
   fChain->SetBranchAddress("gamma_size", &gamma_size, &b_gamma_size);
   fChain->SetBranchAddress("gamma_pt", gamma_pt, &b_gamma_pt);
   fChain->SetBranchAddress("gamma_eta", gamma_eta, &b_gamma_eta);
   fChain->SetBranchAddress("gamma_phi", gamma_phi, &b_gamma_phi);
   fChain->SetBranchAddress("gamma_mass", gamma_mass, &b_gamma_mass);
   fChain->SetBranchAddress("gamma_idvar", gamma_idvar, &b_gamma_idvar);
   fChain->SetBranchAddress("gamma_reliso", gamma_reliso, &b_gamma_reliso);
   fChain->SetBranchAddress("gamma_idpass", gamma_idpass, &b_gamma_idpass);
   fChain->SetBranchAddress("gamma_isopass", gamma_isopass, &b_gamma_isopass);
   fChain->SetBranchAddress("elec_size", &elec_size, &b_elec_size);
   fChain->SetBranchAddress("elec_pt", elec_pt, &b_elec_pt);
   fChain->SetBranchAddress("elec_eta", elec_eta, &b_elec_eta);
   fChain->SetBranchAddress("elec_phi", elec_phi, &b_elec_phi);
   fChain->SetBranchAddress("elec_mass", elec_mass, &b_elec_mass);
   fChain->SetBranchAddress("elec_charge", elec_charge, &b_elec_charge);
   fChain->SetBranchAddress("elec_idvar", elec_idvar, &b_elec_idvar);
   fChain->SetBranchAddress("elec_reliso", elec_reliso, &b_elec_reliso);
   fChain->SetBranchAddress("elec_idpass", elec_idpass, &b_elec_idpass);
   fChain->SetBranchAddress("elec_isopass", elec_isopass, &b_elec_isopass);
   fChain->SetBranchAddress("muon_size", &muon_size, &b_muon_size);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_mass", muon_mass, &b_muon_mass);
   fChain->SetBranchAddress("muon_charge", muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_idvar", muon_idvar, &b_muon_idvar);
   fChain->SetBranchAddress("muon_reliso", muon_reliso, &b_muon_reliso);
   fChain->SetBranchAddress("muon_idpass", muon_idpass, &b_muon_idpass);
   fChain->SetBranchAddress("muon_isopass", muon_isopass, &b_muon_isopass);
   fChain->SetBranchAddress("tau_size", &tau_size, &b_tau_size);
   fChain->SetBranchAddress("tau_pt", tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("tau_eta", tau_eta, &b_tau_eta);
   fChain->SetBranchAddress("tau_phi", tau_phi, &b_tau_phi);
   fChain->SetBranchAddress("tau_mass", tau_mass, &b_tau_mass);
   fChain->SetBranchAddress("tau_charge", tau_charge, &b_tau_charge);
   fChain->SetBranchAddress("tau_decaymode", tau_decaymode, &b_tau_decaymode);
   fChain->SetBranchAddress("tau_neutraliso", tau_neutraliso, &b_tau_neutraliso);
   fChain->SetBranchAddress("tau_chargediso", tau_chargediso, &b_tau_chargediso);
   fChain->SetBranchAddress("tau_combinediso", tau_combinediso, &b_tau_combinediso);
   fChain->SetBranchAddress("tau_isopass", tau_isopass, &b_tau_isopass);
   fChain->SetBranchAddress("jetpuppi_size", &jetpuppi_size, &b_jetpuppi_size);
   fChain->SetBranchAddress("jetpuppi_pt", jetpuppi_pt, &b_jetpuppi_pt);
   fChain->SetBranchAddress("jetpuppi_eta", jetpuppi_eta, &b_jetpuppi_eta);
   fChain->SetBranchAddress("jetpuppi_phi", jetpuppi_phi, &b_jetpuppi_phi);
   fChain->SetBranchAddress("jetpuppi_mass", jetpuppi_mass, &b_jetpuppi_mass);
   fChain->SetBranchAddress("jetpuppi_idpass", jetpuppi_idpass, &b_jetpuppi_idpass);
   fChain->SetBranchAddress("jetpuppi_DeepJET", jetpuppi_DeepJET, &b_jetpuppi_DeepJET);
   fChain->SetBranchAddress("jetpuppi_btag", jetpuppi_btag, &b_jetpuppi_btag);
   fChain->SetBranchAddress("jetchs_size", &jetchs_size, &b_jetchs_size);
   fChain->SetBranchAddress("jetchs_pt", &jetchs_pt, &b_jetchs_pt);
   fChain->SetBranchAddress("jetchs_eta", &jetchs_eta, &b_jetchs_eta);
   fChain->SetBranchAddress("jetchs_phi", &jetchs_phi, &b_jetchs_phi);
   fChain->SetBranchAddress("jetchs_mass", &jetchs_mass, &b_jetchs_mass);
   fChain->SetBranchAddress("jetchs_idpass", &jetchs_idpass, &b_jetchs_idpass);
   fChain->SetBranchAddress("jetchs_DeepJET", &jetchs_DeepJET, &b_jetchs_DeepJET);
   fChain->SetBranchAddress("jetchs_btag", &jetchs_btag, &b_jetchs_btag);
   fChain->SetBranchAddress("fatjet_size", &fatjet_size, &b_fatjet_size);
   fChain->SetBranchAddress("fatjet_pt", fatjet_pt, &b_fatjet_pt);
   fChain->SetBranchAddress("fatjet_eta", fatjet_eta, &b_fatjet_eta);
   fChain->SetBranchAddress("fatjet_phi", fatjet_phi, &b_fatjet_phi);
   fChain->SetBranchAddress("fatjet_mass", fatjet_mass, &b_fatjet_mass);
//    fChain->SetBranchAddress("fatjet_pt", fatjet_pt, &b_fatjet_pt);
   fChain->SetBranchAddress("fatjet_tau1", fatjet_tau1, &b_fatjet_tau1);
   fChain->SetBranchAddress("fatjet_tau2", fatjet_tau2, &b_fatjet_tau2);
   fChain->SetBranchAddress("fatjet_tau3", fatjet_tau3, &b_fatjet_tau3);
   fChain->SetBranchAddress("fatjet_tau4", fatjet_tau4, &b_fatjet_tau4);
   fChain->SetBranchAddress("fatjet_msoftdrop", fatjet_msoftdrop, &b_fatjet_msoftdrop);
   fChain->SetBranchAddress("metpuppi_size", &metpuppi_size, &b_metpuppi_size);
   fChain->SetBranchAddress("metpuppi_pt", metpuppi_pt, &b_metpuppi_pt);
   fChain->SetBranchAddress("metpuppi_phi", metpuppi_phi, &b_metpuppi_phi);
   fChain->SetBranchAddress("metpf_size", &metpf_size, &b_metpf_size);
   fChain->SetBranchAddress("metpf_pt", &metpf_pt, &b_metpf_pt);
   fChain->SetBranchAddress("metpf_phi", &metpf_phi, &b_metpf_phi);
   Notify();
}

Bool_t TauClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TauClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TauClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TauClass_cxx
