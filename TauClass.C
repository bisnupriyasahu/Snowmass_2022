#define TauClass_cxx
#include "TauClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <vector>
#include "TStopwatch.h"
#include <cstring>
#include <list>
#include <TStyle.h>
#include <math.h>

#include "TVector3.h"
using namespace std;


int main(int argc, const char* argv[])
{
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
    {
      std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
      return 1;
    }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
    {
      std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
      return 1;
    }
  TauClass t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;

}
void TauClass::Loop(Long64_t maxEvents, int reportEvery)
{

   if (fChain == 0) return;

   int nTotal;
   nTotal = 0;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "Total entries: " << nentries << std::endl;
   Long64_t nentriesToCheck = nentries;

   if (maxEvents != -1LL && nentries > maxEvents)
     nentriesToCheck = maxEvents;
   nTotal = nentriesToCheck;

   Long64_t nbytes = 0, nb = 0;
   std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
   TStopwatch sw;
   sw.Start();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;                                                                                                                                                                    
     //     for( int i=0; i<nMC; i++ )
     //{}
     int  ngoodtaus = 0;

     int index_tau1, index_tau2;
     index_tau1 = index_tau1 = -99;
     
     std::vector<int> index_btagjet;
     index_btagjet.clear();
     
     index_tau1 = TauCand1();
     index_tau2 = TauCand2(index_tau1);
     
     index_btagjet = BTagjet();
     double ht=0;
     //      std::cout<<"index_btagjet.size()  : "<<index_btagjet.size()<<std::endl;
     //     std::cout<<"met_size is : "<<metpf_size<<std::endl;
     for(int i=0; i<jetpuppi_size; i++)
       {
	 //	 std::cout<<"cominginside 1"<<std::endl;
	 if(jetpuppi_pt[i] > 30 )
	   {
	     //  std::cout<<"cominginside ht loop"<<std::endl;
	     ht+= jetpuppi_pt[i];
	   }
	 
	 for(int i=0; i<tau_size; i++)
	   {
	     
	     //	 if(tau_pt[i] > 30 && fabs(tau_eta)[i] < 3 && bool( jet.tau_isopass() & (1 << 2) ) == True)
	     //if(tau_pt[i] > 30 && fabs(tau_eta[i]) < 3 && bool(tau_isopass[i] & (1 << 2)) == 1)
	     
	     if( index_tau1 > -1 && index_tau2 > -1 )
	       {
		 
		 //std::cout<<"charge of index_tau1: "<<tau_charge[index_tau1]<<std::endl;
		 //std::cout<<"charge of index_tau2: "<<tau_charge[index_tau2]<<std::endl;
		 if(tau_charge[index_tau1]*tau_charge[index_tau2] < 0  )
		   {
		     
		     ngoodtaus++;
		     //std::cout<<"cominginside"<<std::endl;
		
		     
		     // std::cout<<"met_size is : "<<metpf_size<<std::endl;
		     if(index_btagjet.size()>0 && ht > 50)
		       {
			 if(metpuppi_pt[i] > 50)
			   { 
			     tau_pt_->Fill(tau_pt [i]);
			     tau_eta_->Fill(tau_eta [i]);
			     met_pt_->Fill(metpuppi_pt[i]);
			     met_phi_->Fill(metpuppi_phi[i]);


			     std::cout<<"met_pt is : "<<metpuppi_pt[i]<<std::endl;
			     //std::cout<<"coming inside met loop"<<std::endl;
			     
			   }
		       }
		     
		   }
	       }
	   }
       }//tau ptand etacut loop



















	 //	 float min_dR = 999.9;
	 /* for(int j=0; j<genpart_size; j++)
	   {
	     //	     	     std::cout<<"genpart_status [j]"<<genpart_status [j]<<std::endl;
	     float gen_Id;
	     gen_Id = genpart_pid [j];	   
	     if(fabs(gen_Id) == 15)
	       {//	     if(fabs(genpart_pid [j]) = 15 && genpart_status [j] = 1)
	       
		 std::cout<<"coming inside gen"<<std::endl;
		 float dR = sqrt(  ((tau_eta [i] - genpart_eta [j])*(tau_eta [i] - genpart_eta [j])) + ((tau_phi [i] - genpart_phi [j])*(tau_phi [i] - genpart_phi [j]) ) );
		 //float dR = DeltaR(tau_eta [i],tau_phi [i],genpart_eta [j],genpart_phi [j]); 
		 if (dR<min_dR)
		   {
		     min_dR=dR;
		   }
	       }//if condition on taus
	   }//loop of j on gen taus
	 if (min_dR <= 0.15) 
	   {
	     std::cout<<"coming inside dr"<<std::endl;
	     ntau++;
	     tau_pt_->Fill(tau_pt [i]);
	     tau_eta_->Fill(tau_eta [i]);
	 
	     std::cout<<"tau pt is :  "<<tau_pt [i]<<std::endl;	 
	   }
	 //	 std::cout<<"ntau is :  "<<ntau<<std::endl;	 
	 */
	 // std::cout<<"genpart_size  :  "<<genpart_size<<std::endl;
     //}// loop over i Reco taus
     

   }//Jentry
   if((nentriesToCheck-1)%reportEvery != 0)
     std::cout<<"Finished entry "<<(nentriesToCheck-1)<<"/"<<(nentriesToCheck-1)<<std::endl;
   sw.Stop();



}//void Loop

int TauClass::TauCand1()
{

  std::vector<int> tmpCand;
  tmpCand.clear();
  TVector3 tmp_tau;
  for(int i=0; i<tau_size; i++)
    {
      tmp_tau.SetPtEtaPhi(tau_pt[i],tau_eta[i],tau_phi[i]);		      
      if(tau_pt[i] > 30 && fabs(tau_eta[i]) < 3 && tau_isopass[i] >> 2 & 1==1)
	{
	  if(fabs(tau_charge[i]) == 1)
	    {	  tmpCand.push_back(i); }
	}
    }

  if(tmpCand.size()>0)
    return tmpCand[0];
  else
    return -1;
}


std::vector<int> TauClass::BTagjet()
{

  std::vector<int> tmpCand;
  tmpCand.clear();
  for(int i=0; i<jetpuppi_size; i++) 

    {
	  //	  std::cout<<"coming in jet eta: "<<<<std::endl;
	  if(jetchs_pt[i] > 30 && fabs(jetchs_eta[i]) < 5 && jetpuppi_btag[i]>> 1 & 1==1)
	{
    	  //std::cout<<"coming in jet pt: "<<jetpuppi_pt[i]<<std::endl;
	  //	  std::cout<<"coming in jet eta: "<<jetpuppi_eta[i]<<std::endl;
	  tmpCand.push_back(i);

	}
    }

  return tmpCand;

}


int TauClass::TauCand2(int tau1Index)
{
  std::vector<int> tmpCand;
  tmpCand.clear();
  TVector3 tmp_tau_2;
  for(int i=0; i<tau_size; i++)
    {
      if(tau1Index >= 0 && i == tau1Index) continue;
      tmp_tau_2.SetPtEtaPhi(tau_pt[i],tau_eta[i],tau_phi[i]);
      if(tau_pt[i] > 30 && fabs(tau_eta[i]) < 3 && tau_isopass[i] >> 2 & 1==1)
	{if(fabs(tau_charge[i]) == 1)
            {     tmpCand.push_back(i); }
        }
    }
  if(tmpCand.size()>0)
    return tmpCand[0];
  else
    return -1;
}




void TauClass::Histos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  fileName->cd();
  tau_pt_ = new TH1F("tau_pt_", "pt", 15, 0, 5000); 
  tau_eta_ = new TH1F("tau_eta_","eta",20,-5,5);
  met_pt_ = new TH1F("met_pt_", "pt", 15, 0, 5000); 
  met_phi_ = new TH1F("met_phi_","phi",20,-5,5);
  tauVis_pt_ = new TH1F("tauVis_pt_", "pt", 15, 0, 5000); 
  tauVis_eta_ = new TH1F("tauVis_eta_","eta",20,-5,5);
  tauVis_M_ = new TH1F("tauVis_M_", "Mass", 50, 0, 10);
  pi_pt_ = new TH1F("pi_pt_", "pt", 15, 0, 5000);
  ptratio_ = new TH1F("ptratio", "pt", 50, 0, 1.10); 
  tauvispt_taupt_ = new TH1F("tauvispt_taupt_", "pt", 50, 0, 20); 
  pipt_tauvispt_ = new TH1F("pipt_tauvispt_", "pt", 50, 0, 2);

}

