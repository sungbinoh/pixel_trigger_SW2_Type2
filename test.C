#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TFile.h>
#include <TChain.h>

using namespace std;

map<TString, TH1*> maphist;

/// Making Histogram Function //////////////////////////////////////////////// 
void MakeHistograms(TString hname, int nbins, float xmin, float xmax){

  maphist[hname] = new TH1F(hname.Data(), hname.Data(), nbins, xmin, xmax);

}
//////////////////////////////////////////////////////////////////////////////



/// Getting Histogram Function /////////////////////////////////////////////// 
TH1 * GetHist(TString hname){

  TH1 *h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;

  return h;

}
//////////////////////////////////////////////////////////////////////////////



/// Filling Histogram Function ///////////////////////////////////////////////
void FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins){

  if(GetHist(histname)) GetHist(histname) -> Fill(value, w);
  else{
    cout << "Making histogram..." <<endl;
    MakeHistograms(histname, nbins, xmin, xmax);
    if(GetHist(histname)) GetHist(histname) -> Fill(value, w);
  }

}
//////////////////////////////////////////////////////////////////////////////



/// Calculating delta Phi Function ///////////////////////////////////////////    
float delta_phi( float phi1, float phi2){
  float sub_phi = phi2 - phi1;
  float del_phi;
  if( cos(sub_phi / 2) >= 0){
    return sub_phi;
  }
  else if( sub_phi < 0){
    del_phi = sub_phi + 2 * TMath::Pi();
    return del_phi;
  }
  else{
    del_phi = sub_phi - 2 * TMath::Pi();
    return del_phi;
  }
}
//////////////////////////////////////////////////////////////////////////////




void test::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L test.C
//      Root > test t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   //float EgEt, EgEta, EgPhi, EgN;
   //float GenElPT, GenElPhi;
   
   TH2F *EtEt = new TH2F("genEt & EgEt", "GenEt & EgEt", 80, 0, 80, 80, 0, 80);

   //TH2F *GenEt_Dphi = new TH2F("GenEt_Dphi","GenEt vs Dphi", 150, 0., 150., 500, -0.3, 0.3);
   TH2F *GenEt_Dphi = new TH2F("GenEt_Dphi","GenEt vs Dphi", 500, -0.3, 0.3, 60, 0., 60.);
   
   Long64_t nentries = fChain->GetEntriesFast();
   
   //vectors to use to draw histograms
   TVector3 L1_hit;
   TVector3 L2_hit;
   TVector3 L3_hit;
   TVector3 L4_hit;
   TVector3 D1_hit;
   TVector3 D2_hit;
   TVector3 D3_hit;

   int L1L2 = 0;
   int L1L3 = 0;
   int L1L4 = 0;
   int L2L3 = 0;
   int L2L4 = 0;
   int L3L4 = 0;

   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry < nentries ; jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     
     //Vectors to pushback
     std::vector<TVector3> L1_hits;
     std::vector<TVector3> L2_hits;
     std::vector<TVector3> L3_hits;
     std::vector<TVector3> L4_hits;
     std::vector<TVector3> D1_hits;
     std::vector<TVector3> D2_hits;
     std::vector<TVector3> D3_hits;
     std::vector<int> hitted_layers;
     

     float EgN = egCrysN;
     float EgEt, EgEta, EgPhi, EgGx, EgGy, EgGz;
     float GenElPT, GenElPhi, GenElEta;
     int junk;
     //hit number for each pixel layer
     int L1_N = 0;
     int L2_N = 0;
     int L3_N = 0;
     int L4_N = 0;
     int D1_N = 0;
     int D2_N = 0;
     int D3_N = 0;
     
     GenElPT = genPartPt->at(0);
     GenElPhi = genPartPhi->at(0);
     GenElEta = genPartEta->at(0);

     //if(GenElPT > 140) continue;
     
     vector<float> EGET;
     
     float mini = 999.;
     for(int k = 0; k < EgN; k++){
       EgEt = egCrysEt -> at(k);
       if(1){
	 EGET.push_back(EgEt);
	 sort(EGET.begin(), EGET.end(), greater<float>());
	 if(mini >= fabs(GenElPT - EGET[k])){
	   mini = fabs(GenElPT - EGET[k]);
	 }
       }
     }


     int egN_10 = 0;
     for(int q = 0; q < egCrysN; q++){
       if(egCrysEt -> at(q) > 10) egN_10 ++;
     }

     
     if(egN_10 != 1) continue;
     
     for(int q = 0; q < egCrysN; q++){
       
       EgEt = egCrysEt -> at(q);
       EgEta = egCrysEta -> at(q);
       EgPhi = egCrysPhi -> at(q);
       EgGx = egCrysGx -> at(q);
       EgGy = egCrysGy -> at(q);
       EgGz = egCrysGz -> at(q);

       TVector3 eg_vector;
       eg_vector.SetXYZ( EgGx, EgGy, EgGz);
       
       //EgEt cut
       if(EgEt < 10) continue;
       /*
       if(EgEt>= 9 && EgEt<11 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.4)  /(-0.02-0.02)  *(GenElPhi-EgPhi +0.02  )  -0.2 ) continue; }
       else if(EgEt>=11 && EgEt<15 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.4)  /(-0.02-0.02)  *(GenElPhi-EgPhi +0.02  )  -0.2 ) continue; }
       else if(EgEt>=15 && EgEt<20 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.17) / -0.011       *(GenElPhi-EgPhi        )  -0.2 ) continue; }
       else if(EgEt>=20 && EgEt<30 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.17) /(0.018-0.026) *(GenElPhi-EgPhi -0.018 )  -0.2 ) continue; }
       else if(EgEt>=30 && EgEt<40 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.06) /(0.015-0.036) *(GenElPhi-EgPhi -0.015 )  -0.2 ) continue; }
       else if(EgEt>=40 && EgEt<50 ) { if( ((GenElPT-EGET[q])/GenElPT) < (-0.2-0.03) /(0.020-0.042) *(GenElPhi-EgPhi -0.020 )  -0.2 ) continue; }
       else if(EgEt>=50 ) { if( ((GenElPT-EGET[q])/GenElPT) < -0.14 ) continue; }       
       */

       float ratio = EgEt / GenElPT ;
       //if( fabs(ratio - 1) > 0.5 ) continue;
       
       EtEt -> Fill( GenElPT, EgEt );

       float match_R = sqrt( (EgEta - GenElEta) * (EgEta - GenElEta) + (EgPhi - GenElPhi) * (EgPhi -GenElPhi) );

       if(match_R > 0.3) continue;
 
       //to use only R4 region
       //if( fabs(EgEta) < 1.9 || fabs(EgEta) > 2.5 ) continue;
       
       
       for(int bhit = 0; bhit < bHitN; bhit ++){
	 
	 float b_cl_r = sqrt( pow(bHitGx->at(bhit), 2) + pow(bHitGy->at(bhit), 2) );
	 
	 //if(bHitGz->at(bhit) > 28) continue; //z-length cut
	 
	 if(b_cl_r > 2 && b_cl_r < 5 && fabs(bHitGz->at(bhit)) < 28){ 
	   L1_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
	   L1_N ++;
	 }
	 else if(b_cl_r > 6 && b_cl_r < 8 && fabs(bHitGz->at(bhit)) < 28){
	   L2_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
           L2_N ++;
	 }
	 else if(b_cl_r> 10 && b_cl_r < 12 && fabs(bHitGz->at(bhit)) < 28){
	   L3_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
           L3_N ++;
	 }
	 else if(b_cl_r> 14 && b_cl_r < 18 && fabs(bHitGz->at(bhit)) < 28){
	   L4_hits.push_back( TVector3(bHitGx->at(bhit), bHitGy->at(bhit), bHitGz->at(bhit)));
           L4_N ++;
	 }

	 
       }//bhit for
       
       for(int fhit = 0; fhit < fHitN; fhit ++){
	 
	 float f_z = fHitGz->at(fhit);
	 
	 if( fabs(f_z) > 28 && fabs(f_z) < 36){
	   D1_hits.push_back( TVector3(fHitGx->at(fhit), fHitGy->at(fhit), fHitGz->at(fhit)));
	   D1_N ++;
	 }
	 else if( fabs(f_z) > 36 && fabs(f_z) < 44){
	   D2_hits.push_back( TVector3(fHitGx->at(fhit), fHitGy->at(fhit), fHitGz->at(fhit)));
	   D2_N ++;
	 }
	 else if( fabs(f_z) > 44 && fabs(f_z) < 54){
	   D3_hits.push_back( TVector3(fHitGx->at(fhit), fHitGy->at(fhit), fHitGz->at(fhit)));
	   D3_N ++;
	 }
	 else junk = 1;
	 
       }//fhit for
       
       
       for(int range = 10; range < 141; range ++){ 
	 if(EgEt >= range && EgEt < range+1 ){
	   char range_in_char[30];//TString range_in_string;   
	   int range_sub = range;
	   sprintf( range_in_char, "%d", range_sub);
	   TString range_in_string = range_in_char;
	   if(range <10){
	     range_in_string.Insert(0, "00");
	   }
	   if(range >9 && range < 100){
	     range_in_string.Insert(0, "0");
	   }
	   TString histname_phi = "_layer_delta_phi_";
	   TString histname_eta = "_layer_delta_eta_";
	   TString histname_R = "_layer_delta_R_";
	   
	   histname_phi.Append(range_in_string);
	   histname_eta.Append(range_in_string);
	   histname_R.Append(range_in_string);

	   TString L1L2_string = "L1L2_";
	   TString L1L3_string = "L1L3_";
	   TString L1L4_string = "L1L4_";
	   TString L2L3_string = "L2L3_";
	   TString L2L4_string = "L2L4_";
	   TString L3L4_string = "L3L4_";

	   TString L1D1_string = "L1D1_";
	   TString L2D1_string = "L2D1_";
	   TString L3D1_string = "L3D1_";

	   TString L1D2_string = "L1D2_";
	   TString L2D2_string = "L2D2_";
	   TString D1D2_string = "D1D2_";

	   TString L1D3_string = "L1D3_";
	   TString D1D3_string = "D1D3_";
	   TString D2D3_string = "D2D3_";
	   
	   /*
	   //Region 1, L1234
	   
	   if( fabs(EgEta) < 1.3 ){
	    
	     L1L2_string.Insert(0, "R1_");
	     L1L3_string.Insert(0, "R1_");
	     L1L4_string.Insert(0, "R1_");
	     L2L3_string.Insert(0, "R1_");
	     L2L4_string.Insert(0, "R1_");
	     L3L4_string.Insert(0, "R1_");

	     if(L1_N > 0 && L2_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
		 for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){

		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];

		   float current_dPhi = delta_phi( current_L2.Phi(), current_L1.Phi() );
		   float current_dEta = current_L2.Eta() - current_L1.Eta();
		   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
		   FillHist(L1L2_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
		   FillHist(L1L2_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
		   FillHist(L1L2_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
		 }
	       }
	     }
	     if(L1_N > 0 && L3_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
		 for(int hitnum_L3 = 0; hitnum_L3 < L3_N; hitnum_L3++){

		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_L3 = eg_vector - L3_hits[hitnum_L3];

                   float current_dPhi = delta_phi( current_L3.Phi(), current_L1.Phi() );
                   float current_dEta =current_L3.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

		   FillHist(L1L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
		   FillHist(L1L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
		   FillHist(L1L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
		 }
	       }
	     }
	     if(L1_N > 0 && L4_N > 0){
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
		 for(int hitnum_L4 = 0; hitnum_L4 < L4_N; hitnum_L4++){

		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_L4 = eg_vector - L4_hits[hitnum_L4];

                   float current_dPhi = delta_phi( current_L4.Phi(), current_L1.Phi() );
                   float current_dEta =current_L4.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

		   FillHist(L1L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
		   FillHist(L1L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
		   FillHist(L1L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
		 }
	       }
	     }
	     if(L2_N > 0 && L3_N > 0){
	       for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
		 for(int hitnum_L3 = 0; hitnum_L3 < L3_N; hitnum_L3++){

		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];
                   TVector3 current_L3 = eg_vector - L3_hits[hitnum_L3];

                   float current_dPhi = delta_phi( current_L3.Phi(), current_L2.Phi() );
                   float current_dEta =current_L3.Eta() - current_L2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

		   FillHist(L2L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
		   FillHist(L2L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
		   FillHist(L2L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
		 }
	       }
	     }
	     if(L2_N > 0 && L4_N > 0){
	       for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
		 for(int hitnum_L4 = 0; hitnum_L4 < L4_N; hitnum_L4++){
		   
		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];
                   TVector3 current_L4 = eg_vector - L4_hits[hitnum_L4];

                   float current_dPhi = delta_phi( current_L4.Phi(), current_L2.Phi() );
                   float current_dEta =current_L4.Eta() - current_L2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
		   FillHist(L2L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
		   FillHist(L2L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
		   FillHist(L2L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
		 }
	       }
	     }
	     if(L3_N > 0 && L4_N > 0){
	       for(int hitnum_L3 = 0; hitnum_L3 < L3_N; hitnum_L3++){
		 for(int hitnum_L4 = 0; hitnum_L4 < L4_N; hitnum_L4++){
		  
		   TVector3 current_L3 = eg_vector - L3_hits[hitnum_L3];
                   TVector3 current_L4 = eg_vector - L4_hits[hitnum_L4];

                   float current_dPhi = delta_phi( current_L4.Phi(), current_L3.Phi() );
                   float current_dEta =current_L4.Eta() - current_L3.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

		   FillHist(L3L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
		   FillHist(L3L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
		   FillHist(L3L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
		 }
	       }
	     }
	     
	   }//*/
	   
	   
	   /*
	   //Region 2, L123D1
	   if( fabs(EgEta) > 1.3 && fabs(EgEta) < 1.6 ){
	     
	     L1L2_string.Insert(0, "R2_");
             L1L3_string.Insert(0, "R2_");
             L1L4_string.Insert(0, "R2_");
             L2L3_string.Insert(0, "R2_");
             L2L4_string.Insert(0, "R2_");
             L3L4_string.Insert(0, "R2_");

	     
	     if(L1_N > 0 && L2_N > 0){
               L1L2 ++;
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
		   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];

                   float current_dPhi = delta_phi( current_L2.Phi(), current_L1.Phi() );
                   float current_dEta = current_L2.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L1L2_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(L1_N > 0 && L3_N > 0){
               L1L3 ++;
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_L3 = 0; hitnum_L3 < L3_N; hitnum_L3++){
		   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_L3 = eg_vector - L3_hits[hitnum_L3];

                   float current_dPhi = delta_phi( current_L3.Phi(), current_L1.Phi() );
                   float current_dEta = current_L3.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L1L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(L1_N > 0 && D1_N > 0){
               L1L4 ++;
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L1.Phi() );
                   float current_dEta = current_D1.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L1L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(L2_N > 0 && L3_N > 0){
	       L2L3 ++;
	       for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
                 for(int hitnum_L3 = 0; hitnum_L3 < L3_N; hitnum_L3++){
		   
		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];
                   TVector3 current_L3 = eg_vector - L3_hits[hitnum_L3];

                   float current_dPhi = delta_phi( current_L3.Phi(), current_L2.Phi() );
                   float current_dEta = current_L3.Eta() - current_L2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(L2_N > 0 && D1_N > 0){
	       L2L4 ++;
	       for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
		   
		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];
                   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L2.Phi() );
                   float current_dEta = current_D1.Eta() - current_L2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L2L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(L3_N > 0 && D1_N > 0){
	       L3L4 ++;
	       for(int hitnum_L3 = 0; hitnum_L3 < L3_N; hitnum_L3++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){

		   TVector3 current_L3 = eg_vector - L3_hits[hitnum_L3];
                   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L3.Phi() );
                   float current_dEta = current_D1.Eta() - current_L3.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L3L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }

	     
	   }//*/
	   

	   /*
	   //Region 3, L12D12
	   if( fabs(EgEta) > 1.6 && fabs(EgEta) < 1.9 ){

	     L1L2_string.Insert(0, "R3_");
             L1L3_string.Insert(0, "R3_");
             L1L4_string.Insert(0, "R3_");
             L2L3_string.Insert(0, "R3_");
             L2L4_string.Insert(0, "R3_");
             L3L4_string.Insert(0, "R3_");


             if(L1_N > 0 && L2_N > 0){
               L1L2 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){

		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];

                   float current_dPhi = delta_phi( current_L2.Phi(), current_L1.Phi() );
                   float current_dEta = current_L2.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L1L2_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L1_N > 0 && D1_N > 0){
               L1L3 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){

		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L1.Phi() );
                   float current_dEta = current_D1.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L1L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L1_N > 0 && D2_N > 0){
	       L1L4 ++;
	       for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
		   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_L1.Phi() );
                   float current_dEta = current_D2.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L1L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(L2_N > 0 && D1_N > 0){
               L2L3 ++;
               for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
		   
		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];
                   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L2.Phi() );
                   float current_dEta = current_D1.Eta() - current_L2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L2_N > 0 && D2_N > 0){
               L2L4 ++;
               for(int hitnum_L2 = 0; hitnum_L2 < L2_N; hitnum_L2++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
                   
		   TVector3 current_L2 = eg_vector - L2_hits[hitnum_L2];
                   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_L2.Phi() );
                   float current_dEta = current_D2.Eta() - current_L2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(D1_N > 0 && D2_N > 0){
               L3L4 ++;
               for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){

		   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];
                   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_D1.Phi() );
                   float current_dEta = current_D2.Eta() - current_D1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L3L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	   }//*/
	   
	   
	   /*
	   //Region 4, L1D123
	   if( fabs(EgEta) > 1.9 && fabs(EgEta) < 2.5 ){

	     L1L2_string.Insert(0, "R4_");
             L1L3_string.Insert(0, "R4_");
             L1L4_string.Insert(0, "R4_");
             L2L3_string.Insert(0, "R4_");
             L2L4_string.Insert(0, "R4_");
             L3L4_string.Insert(0, "R4_");

	     
             if(L1_N > 0 && D1_N > 0){
               L1L2 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
		   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L1.Phi() );
                   float current_dEta = current_D1.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L1L2_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L1_N > 0 && D2_N > 0){
               L1L3 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
                   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
		   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_L1.Phi() );
                   float current_dEta = current_D2.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L1L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L1_N > 0 && D3_N > 0){
               L1L4 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D3 = 0; hitnum_D3 < D3_N; hitnum_D3++){
                   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D3 = eg_vector - D3_hits[hitnum_D3];

                   float current_dPhi = delta_phi( current_D3.Phi(), current_L1.Phi() );
                   float current_dEta = current_D3.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L1L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(D1_N > 0 && D2_N > 0){
               L2L3 ++;
               for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
                   
		   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];
                   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_D1.Phi() );
                   float current_dEta = current_D2.Eta() - current_D1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(D1_N > 0 && D3_N > 0){
               L2L4 ++;
               for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                 for(int hitnum_D3 = 0; hitnum_D3 < D3_N; hitnum_D3++){
                   
		   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];
                   TVector3 current_D3 = eg_vector - D3_hits[hitnum_D3];

                   float current_dPhi = delta_phi( current_D3.Phi(), current_D1.Phi() );
                   float current_dEta = current_D3.Eta() - current_D1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(D2_N > 0 && D3_N > 0){
               L3L4 ++;
               for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
                 for(int hitnum_D3 = 0; hitnum_D3 < D3_N; hitnum_D3++){
		   
		   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];
                   TVector3 current_D3 = eg_vector - D3_hits[hitnum_D3];

                   float current_dPhi = delta_phi( current_D3.Phi(), current_D2.Phi() );
                   float current_dEta = current_D3.Eta() - current_D2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L3L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }

	   }//*/
	   
	   
	   
	   //Region 5, L1D123
	   if( fabs(EgEta) > 2.5 && fabs(EgEta) < 2.8 ){
	
	     L1L2_string.Insert(0, "R5_");
             L1L3_string.Insert(0, "R5_");
             L1L4_string.Insert(0, "R5_");
             L2L3_string.Insert(0, "R5_");
             L2L4_string.Insert(0, "R5_");
             L3L4_string.Insert(0, "R5_");


             if(L1_N > 0 && D1_N > 0){
               L1L2 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];

                   float current_dPhi = delta_phi( current_D1.Phi(), current_L1.Phi() );
                   float current_dEta = current_D1.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );

                   FillHist(L1L2_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L2_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L1_N > 0 && D2_N > 0){
               L1L3 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
                   
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_L1.Phi() );
                   float current_dEta = current_D2.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L1L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(L1_N > 0 && D3_N > 0){
               L1L4 ++;
               for(int hitnum_L1 = 0; hitnum_L1 < L1_N; hitnum_L1++){
                 for(int hitnum_D3 = 0; hitnum_D3 < D3_N; hitnum_D3++){
                  
		   TVector3 current_L1 = eg_vector - L1_hits[hitnum_L1];
                   TVector3 current_D3 = eg_vector - D3_hits[hitnum_D3];

                   float current_dPhi = delta_phi( current_D3.Phi(), current_L1.Phi() );
                   float current_dEta = current_D3.Eta() - current_L1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L1L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L1L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(D1_N > 0 && D2_N > 0){
               L2L3 ++;
               for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                 for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){

		   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];
                   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];

                   float current_dPhi = delta_phi( current_D2.Phi(), current_D1.Phi() );
                   float current_dEta = current_D2.Eta() - current_D1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L3_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L3_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
             if(D1_N > 0 && D3_N > 0){
               L2L4 ++;
               for(int hitnum_D1 = 0; hitnum_D1 < D1_N; hitnum_D1++){
                 for(int hitnum_D3 = 0; hitnum_D3 < D3_N; hitnum_D3++){
                   
		   TVector3 current_D1 = eg_vector - D1_hits[hitnum_D1];
                   TVector3 current_D3 = eg_vector - D3_hits[hitnum_D3];

                   float current_dPhi = delta_phi( current_D3.Phi(), current_D1.Phi() );
                   float current_dEta = current_D3.Eta() - current_D1.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L2L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L2L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }
	     if(D2_N > 0 && D3_N > 0){
               L3L4 ++;
               for(int hitnum_D2 = 0; hitnum_D2 < D2_N; hitnum_D2++){
                 for(int hitnum_D3 = 0; hitnum_D3 < D3_N; hitnum_D3++){
		   
		   TVector3 current_D2 = eg_vector - D2_hits[hitnum_D2];
                   TVector3 current_D3 = eg_vector - D3_hits[hitnum_D3];

                   float current_dPhi = delta_phi( current_D3.Phi(), current_D2.Phi() );
                   float current_dEta = current_D3.Eta() - current_D2.Eta();
                   float current_dR = sqrt( pow(current_dPhi, 2) + pow(current_dEta, 2) );
		   
                   FillHist(L3L4_string + histname_phi, current_dPhi, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_eta, current_dEta, 1., -0.03, 0.03, 1000);
                   FillHist(L3L4_string + histname_R,   current_dR  , 1., -0.01, 0.03, 1000);
                 }
               }
             }


	   }
	   //*/
	   
	   
	 }//if EgEt range
	 
	 
       }//for range (Pt)
       
       
     }//EgN for
     
     
     
     
     
     
     
     
     if( jentry % 10000 == 0){
       cout << jentry << " / " << nentries << "  done" << endl;
     }
     
     
     
     // if (Cut(ientry) < 0) continue;
   }//end of nentries loop
   
   TFile *file = new TFile("First_result.root", "recreate");
   for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit ++){
     mapit->second->Write();
   }
   
   cout << "job is finished  "<< endl;
   
   cout << "L1L2 : " << L1L2 << endl;
   cout << "L1L3 : " << L1L3 << endl;
   cout << "L1L4 : " << L1L4 << endl;
   cout << "L2L3 : " << L2L3 << endl;
   cout << "L2L4 : " << L2L4 << endl;
   cout << "L3L4 : " << L3L4 << endl;


 
   TCanvas * c = new TCanvas("", "", 800, 600);
   c->cd();
   GenEt_Dphi -> Draw();
   c->Update();
   c->SaveAs("./GenEt_Dphi.pdf");


  
}
