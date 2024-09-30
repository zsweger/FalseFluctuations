#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "PDGData.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "stdlib.h"
#include "EnergyConfig.h"

using namespace std;

int main(int argc, char** argv)
{
    
    const char* jobname  = argv[1];
    const char* filelist = argv[2];
    const char* failureRateAsChar = argv[3];
    double failureFraction = std::stod(argv[3]);

    TRandom3 *rand = new TRandom3(0);

    ifstream input(filelist);
    string line;
    vector<string> InputList;
    while(input >> line) {
        InputList.push_back(line);
    }
    
    int refmult3[2] = {0,0};
    int npart=0;
    float b;
    
    TFile* file;
    TTree* tree;
    PDGData pdg;
    int   Entries = 0;
    float px[10000];
    float py[10000];
    float pz[10000];
    int   pid[10000];
    int   mul=0;
    
    int np[29];

 
    TFile *output = new TFile(Form("%s_%sfailureRateAsChar.root",jobname,failureRateAsChar),"RECREATE");
    TTree *otree = new TTree("tree","UrQMD Event Tree");
    otree->Branch("Npart",   &npart,    "Npart/I");
    otree->Branch("b",       &b,        "b/f");
    otree->Branch("refmult3",&refmult3, "refmult3[2]/I");
    //np[0]: Standard no pileup
    //np[1]: Standard pileup with pileup protons and pileup FXTMult3
    //np[2]: Nonstandard pileup with no pileup protons and pileup FXTMult3
    //np[3]: Nonstandard pileup with half pileup protons and pileup FXTMult3
    otree->Branch("np",       np,       "np[29]/I");// first loop for pT sencond lop for rapidity
    
    TH2D *hpT_y = new TH2D("hpT_y", ";Rapidity;p_{T}[GeV/c]", 400, -2, 2, 500, 0, 5);
    TH1D *refMult3[2];
    TH1F *impactParameter = new TH1F("impactParameter","impact parameter",200,0,20);
    TH1F *npBB = new TH1F("npBB","npBB",100,-0.5,99.5);
    TH1F *npBA = new TH1F("npBA","npBA",100,-0.5,99.5);
    TH1F *npAB = new TH1F("npAB","npAB",100,-0.5,99.5);
    TH1F *npAA = new TH1F("npAA","npAA",100,-0.5,99.5);
    TH1F *acceptanceFraction = new TH1F("acceptanceFraction","acceptanceFraction",50,-0.5,49.5);
    TH2F* refMult3AndProtons_all = new TH2F("refMult3AndProtons_all","all",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_allgodgod = new TH2F("refMult3AndProtons_allgodgod","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_allgodbad = new TH2F("refMult3AndProtons_allgodbad","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_allbadgod = new TH2F("refMult3AndProtons_allbadgod","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_allbadbad = new TH2F("refMult3AndProtons_allbadbad","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_godgod = new TH2F("refMult3AndProtons_godgod","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_godbad = new TH2F("refMult3AndProtons_godbad","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_badgod = new TH2F("refMult3AndProtons_badgod","",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_badbad = new TH2F("refMult3AndProtons_badbad","",500,-0.5,499.5,200,-0.5,199.5);
    for(int jm=0;jm<2;++jm){
            refMult3[jm] = new TH1D(Form("refMult3_%d",jm),"refMult3",800,-0.5,799.5);
    }

    cout<<"Starting .."<<endl;
    // root file loop
    for(int i = 0; i<InputList.size(); i++) {
        cout<<i+1<<" file finished .."<<endl;
        file = TFile::Open(InputList[i].c_str());
        file->GetObject("urqmd", tree);
        tree->SetBranchAddress("px",    px);
        tree->SetBranchAddress("py",    py);
        tree->SetBranchAddress("pz",    pz);
        tree->SetBranchAddress("pid",   pid);
        tree->SetBranchAddress("mul",  &mul);
        tree->SetBranchAddress("Npart",&npart);
        tree->SetBranchAddress("b",&b);
        // event loop
        Entries = tree->GetEntries();
        for (int j = 0; j< Entries; j++) {
            tree->GetEntry(j);
            if(npart == 0) continue;
            for(int jm=0;jm<29;++jm){
               np[jm] = 0;
            }

	   //double failureFraction = 0.00001;
	   //double failureFraction = 0.00002;
	   //double failureFraction = 0.00005;
    	   //double failureFraction = 0.0001;
    	   //double failureFraction = 0.0002;
    	   //double failureFraction = 0.0005;
    	   //double failureFraction = 0.001;
           //double failureFraction = 0.002;
           //double failureFraction = 0.005; 
           //double failureFraction = 0.01; //default
	   //double failureFraction = 0.02;
	   //double failureFraction = 0.05;
	   //double failureFraction = 0.10;

	   bool isgood = ((rand->Rndm())>failureFraction);
           float eta_tpc_loedge = -2.15; //approximate upper edge of tpc eta in FXT lab frame
           float eta_tpc_hiedge  = 0.0;
           float y_cm = 0.0;
	   if(ENERGY==3.2) y_cm = -1.139;
	   if(ENERGY==3.5) y_cm = -1.254;
           if(ENERGY==3.9) y_cm = -1.375;
           if(ENERGY==4.5) y_cm = -1.522;
           if(ENERGY==5.2) y_cm = -1.683;
           if(ENERGY==6.2) y_cm = -1.867;
           if(ENERGY==7.2) y_cm = -2.021;
           if(ENERGY==7.7) y_cm = -2.102;
	   if(!isgood)cout<<"Not good!"<<endl;
            refmult3[0] = 0;
            refmult3[1] = 0;
            for(int k=0; k<mul; k++) {
                int apid     = abs(pid[k]);
                int charge   = (pid[k] > 0 ? 1 : -1);
                int bcharge  = (pdg[apid].charge / 3) * (pid[k] > 0 ? 1 : -1);
                float m0     = pdg[apid].m0;
                float p      = sqrt( px[k]*px[k] + py[k]*py[k] + pz[k]*pz[k] );
                float pt     = sqrt( px[k]*px[k] + py[k]*py[k] );
		float phi    = acos(px[k]/pt);
          	if(py[k]<0) phi = -1.0*phi;
                float eta    = 0.5*log( (p+pz[k]) / (p-pz[k]) );
                float aeta   = fabs(eta);
                float E      = sqrt( p*p + m0*m0 );
                float y      = 0.5*log( (E+pz[k]) / (E-pz[k]) );
                float ay     = fabs(y);
		float lab_rap = y_cm-y;
		float lab_p  = 0.5*sqrt(pt*pt*exp(-2.0*lab_rap) + pt*pt*exp(2.0*lab_rap) + 2.0*pt*pt + m0*m0*exp(-2.0*lab_rap) + m0*m0*exp(2.0*lab_rap) - 2.0*m0*m0);
		float y_tpc_loedge = log((sqrt(m0*m0+pt*pt*cosh(eta_tpc_loedge)*cosh(eta_tpc_loedge))+pt*sinh(eta_tpc_loedge))/(sqrt(m0*m0+pt*pt)));
		//float y_tpc_hiedge = log((sqrt(m0*m0+pt*pt*cosh(eta_tpc_hiedge)*cosh(eta_tpc_hiedge))+pt*sinh(eta_tpc_hiedge))/(sqrt(m0*m0+pt*pt)));
		float y_tpc_hiedge = 0.0;//same as above
		bool isgoodsector = true;
		if(phi<0 && !isgood) isgoodsector = false;
		if( (pt > 0.06) && (lab_rap > y_tpc_loedge) && (lab_rap < y_tpc_hiedge) && (apid == 211 || apid == 321)) {
//		if( pt < 2.0  &&  pt > 0.4 && y>-0.5  && y<0.0 && (apid == 211 || apid == 321)) {
			refmult3[0]++;
			if(isgoodsector) refmult3[1]++;
		}
                //==================
                if(pid[k] == 2212){//} && y<0. && y>-0.5 ){
                    
                     if(pt > 2.0  || pt < 0.4 ) continue;
                     if(y  > 0.0  || y  < -0.5) continue;
                     //here plot proton accrptance just to test the acceptance is right
                     hpT_y -> Fill( y, pt);
                     np[0]++;
		     np[2]++;
	 	     if(isgoodsector) {np[1]++;np[3]++;}
		     acceptanceFraction->Fill(40);
                     if(fabs(lab_p)<1.0){ acceptanceFraction->Fill(4);}
                     if(fabs(lab_p)<1.1){ acceptanceFraction->Fill(5);}
                     if(fabs(lab_p)<1.2){ acceptanceFraction->Fill(6);}
                     if(fabs(lab_p)<1.3){ acceptanceFraction->Fill(7);}
                     if(fabs(lab_p)<1.4){ acceptanceFraction->Fill(8);}
                     if(fabs(lab_p)<1.5){ acceptanceFraction->Fill(9);}
                     if(fabs(lab_p)<1.6){ acceptanceFraction->Fill(10);}
                     if(fabs(lab_p)<1.7){ acceptanceFraction->Fill(11);}
                     if(fabs(lab_p)<1.8){ acceptanceFraction->Fill(12);}
                     if(fabs(lab_p)<1.9){ acceptanceFraction->Fill(13);}
                     if(fabs(lab_p)<2.0){ acceptanceFraction->Fill(14);}
                     if(fabs(lab_p)<2.1){ acceptanceFraction->Fill(15);}
                     if(fabs(lab_p)<2.2){ acceptanceFraction->Fill(16);}
                     if(fabs(lab_p)<2.4){ acceptanceFraction->Fill(17);}
                     if(fabs(lab_p)<2.6){ acceptanceFraction->Fill(18);}
                     if(fabs(lab_p)<2.8){ acceptanceFraction->Fill(19);}
                     if(fabs(lab_p)<3.0){ acceptanceFraction->Fill(20);}
                     if(fabs(lab_p)<3.2){ acceptanceFraction->Fill(21);}
                     if(fabs(lab_p)<3.4){ acceptanceFraction->Fill(22);}
                     if(fabs(lab_p)<3.6){ acceptanceFraction->Fill(23);}
                     if(fabs(lab_p)<3.8){ acceptanceFraction->Fill(24);}
                     if(fabs(lab_p)<4.0){ acceptanceFraction->Fill(25);}
                     if(fabs(lab_p)<4.2){ acceptanceFraction->Fill(26);}
                     if(fabs(lab_p)<4.4){ acceptanceFraction->Fill(27);}
                     if(fabs(lab_p)<4.6){ acceptanceFraction->Fill(28);}

                     if(fabs(lab_p)<1.0 || isgoodsector){np[4]++;}
                     if(fabs(lab_p)<1.1 || isgoodsector){np[5]++;}
                     if(fabs(lab_p)<1.2 || isgoodsector){np[6]++;}
                     if(fabs(lab_p)<1.3 || isgoodsector){np[7]++;}
		     if(fabs(lab_p)<1.4 || isgoodsector){np[8]++;}
                     if(fabs(lab_p)<1.5 || isgoodsector){np[9]++;}
                     if(fabs(lab_p)<1.6 || isgoodsector){np[10]++;}
                     if(fabs(lab_p)<1.7 || isgoodsector){np[11]++;}
                     if(fabs(lab_p)<1.8 || isgoodsector){np[12]++;}
                     if(fabs(lab_p)<1.9 || isgoodsector){np[13]++;}
                     if(fabs(lab_p)<2.0 || isgoodsector){np[14]++;}
                     if(fabs(lab_p)<2.1 || isgoodsector){np[15]++;}
                     if(fabs(lab_p)<2.2 || isgoodsector){np[16]++;}
                     if(fabs(lab_p)<2.4 || isgoodsector){np[17]++;}
                     if(fabs(lab_p)<2.6 || isgoodsector){np[18]++;}
                     if(fabs(lab_p)<2.8 || isgoodsector){np[19]++;}
                     if(fabs(lab_p)<3.0 || isgoodsector){np[20]++;}
                     if(fabs(lab_p)<3.2 || isgoodsector){np[21]++;}
                     if(fabs(lab_p)<3.4 || isgoodsector){np[22]++;}
                     if(fabs(lab_p)<3.6 || isgoodsector){np[23]++;}
                     if(fabs(lab_p)<3.8 || isgoodsector){np[24]++;}
                     if(fabs(lab_p)<4.0 || isgoodsector){np[25]++;}
                     if(fabs(lab_p)<4.2 || isgoodsector){np[26]++;}
                     if(fabs(lab_p)<4.4 || isgoodsector){np[27]++;}
                     if(fabs(lab_p)<4.6 || isgoodsector){np[28]++;}

                     if(fabs(lab_p)<2.1 || isgoodsector){np[28]++;}
                }
            } //event loop ends

	    if(refmult3[0]<210 && refmult3[0]>139) npAA->Fill(np[0]);
            if(refmult3[0]<210 && refmult3[0]>139) npAB->Fill(np[1]);
            if(refmult3[1]<210 && refmult3[1]>139) npBA->Fill(np[2]);
            if(refmult3[1]<210 && refmult3[1]>139) npBB->Fill(np[3]);
	    /*if(refmult3[0]<81 && refmult3[0]>56) npAA->Fill(np[0]);
            if(refmult3[0]<81 && refmult3[0]>56) npAB->Fill(np[1]);
            if(refmult3[1]<81 && refmult3[1]>56) npBA->Fill(np[2]);
            if(refmult3[1]<81 && refmult3[1]>56) npBB->Fill(np[3]);*/
            otree->Fill();
	    if(!isgood)cout<<"good,bad proton number: "<<np[2]<<", "<<np[3]<<endl;
            if(!isgood)cout<<"good,bad refmult: "<<refmult3[0]<<", "<<refmult3[1]<<endl;
	    refMult3[0]->Fill(refmult3[0]);
            refMult3[1]->Fill(refmult3[1]);
	    impactParameter->Fill(b);
	    refMult3AndProtons_all->Fill(refmult3[0],np[0]);
	    if(!isgood){
	      refMult3AndProtons_godgod->Fill(refmult3[0],np[0]);
              refMult3AndProtons_godbad->Fill(refmult3[0],np[1]);
              refMult3AndProtons_badgod->Fill(refmult3[1],np[0]);
              refMult3AndProtons_badbad->Fill(refmult3[1],np[1]);
	    }
            refMult3AndProtons_allgodgod->Fill(refmult3[0],np[0]);
            refMult3AndProtons_allgodbad->Fill(refmult3[0],np[1]);
            refMult3AndProtons_allbadgod->Fill(refmult3[1],np[0]);
            refMult3AndProtons_allbadbad->Fill(refmult3[1],np[1]);
        }
        file->Close();
    }// root file loop ends

    output->cd();
    otree->Write();
    TVectorD* v = new TVectorD(1);
    v[0] = failureFraction;
    v->Write("failureFraction");
    hpT_y -> Write();
    acceptanceFraction->Write();
    refMult3[0]->Write();
    refMult3[1]->Write();
    impactParameter->Write();
    npBB->Write();
    npBA->Write();
    npAB->Write();
    npAA->Write();
    refMult3AndProtons_all->Write();
    refMult3AndProtons_allgodgod->Write();
    refMult3AndProtons_allgodbad->Write();
    refMult3AndProtons_allbadgod->Write();
    refMult3AndProtons_allbadbad->Write();
    refMult3AndProtons_godgod->Write();
    refMult3AndProtons_godbad->Write();
    refMult3AndProtons_badgod->Write();
    refMult3AndProtons_badbad->Write();
    output->Close();
    
    return 0;
}

