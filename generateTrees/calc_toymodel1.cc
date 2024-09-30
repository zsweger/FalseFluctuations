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
    const char* pileupRateAsChar = argv[3];
    double pileupFraction = std::stod(argv[3]);


    TRandom3 *rand = new TRandom3(0);
    double pi = 3.14159;

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
    
    int np[27];

    double pileupNp[27]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double pileupRefmult3=0;    

    TFile *output = new TFile(Form("%s_%spileupRate.root",jobname,pileupRateAsChar),"RECREATE");
    TTree *otree = new TTree("tree","UrQMD Event Tree");
    otree->Branch("Npart",   &npart,    "Npart/I");
    otree->Branch("b",       &b,        "b/f");
    otree->Branch("refmult3",&refmult3, "refmult3[2]/I");
    //np[0]: Standard no pileup
    //np[1]: Standard pileup with pileup protons and pileup FXTMult3
    //np[2]: Nonstandard pileup with no pileup protons and pileup FXTMult3
    //np[3]: Nonstandard pileup with half pileup protons and pileup FXTMult3
    otree->Branch("np",       np,       "np[27]/I");

    TH2D *hpT_y = new TH2D("hpT_y", ";Rapidity;p_{T}[GeV/c]", 400, -2, 2, 500, 0, 5);
    TH2D *hpT_y_pileuptest = new TH2D("hpT_y_pileuptest", ";Rapidity;p_{T}[GeV/c]", 400, -2, 2, 500, 0, 5);
    TH2F* refMult3AndProtons_all = new TH2F("refMult3AndProtons_all","all",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_allsinglesingle = new TH2F("refMult3AndProtons_allsinglesingle","all single and single",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_allsingledouble = new TH2F("refMult3AndProtons_allsingledouble","all single and double",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_alldoublesingle = new TH2F("refMult3AndProtons_alldoublesingle","all double and single",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_alldoubledouble = new TH2F("refMult3AndProtons_alldoubledouble","all double and double",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_singlesingle = new TH2F("refMult3AndProtons_singlesingle","single and single",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_singledouble = new TH2F("refMult3AndProtons_singledouble","single and double",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_doublesingle = new TH2F("refMult3AndProtons_doublesingle","double and single",500,-0.5,499.5,200,-0.5,199.5);
    TH2F* refMult3AndProtons_doubledouble = new TH2F("refMult3AndProtons_doubledouble","double and double",500,-0.5,499.5,200,-0.5,199.5);
    TH1D *refMult3[2];
    TH1F *impactParameter = new TH1F("impactParameter","impact parameter",200,0,20);
    TH1F *npBB = new TH1F("npBB","npBB",100,-0.5,99.5);
    TH1F *npBA = new TH1F("npBA","npBA",100,-0.5,99.5);
    TH1F *npAB = new TH1F("npAB","npAB",100,-0.5,99.5);
    TH1F *npAA = new TH1F("npAA","npAA",100,-0.5,99.5);
    TH1F *acceptanceFraction = new TH1F("acceptanceFraction","acceptanceFraction",50,-0.5,49.5);
    for(int jm=0;jm<2;++jm){
            refMult3[jm] = new TH1D(Form("refMult3_%d",jm),"refMult3",800,-0.5,799.5);
    }

    //double pileupFraction = 0.00001;
    //double pileupFraction = 0.00002;
    //double pileupFraction = 0.00005;
    //double pileupFraction = 0.0001;
    //double pileupFraction = 0.0002;
    //double pileupFraction = 0.0005;
    //double pileupFraction = 0.001;
    //double pileupFraction = 0.002; //default
    //double pileupFraction = 0.005; 
    //double pileupFraction = 0.01;

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
	    bool isPileup = ((rand->Rndm())<pileupFraction);
            for(int jm=0;jm<27;++jm){
               np[jm] = 0;
            }
           int np4=0;
           int np5=0;
           int np6=0;
           int np7=0;
           int np8=0;
           int np9=0;
           int np10=0;
           int np11=0;
           int np12=0;
           int np13=0;
           int np14=0;
           int np15=0;
           int np16=0;
           int np17=0;
           int np18=0;
           int np19=0;
           int np20=0;
           int np21=0;
           int np22=0;
           int np23=0;
           int np24=0;
           int np25=0;
           int np26=0;
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
	   //cout<<"y_cm = "<<y_cm<<endl;
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
		phi = phi+pi; // 0 to 2pi
                float eta    = 0.5*log( (6.0*p+pz[k]) / (p-pz[k]) );
                float aeta   = fabs(eta);
                float E      = sqrt( p*p + m0*m0 );
                float y      = 0.5*log( (E+pz[k]) / (E-pz[k]) );
                float ay     = fabs(y);
		float lab_rap = y_cm-y;
		float lab_p  = 0.5*sqrt(pt*pt*exp(-2.0*lab_rap) + pt*pt*exp(2.0*lab_rap) + 2.0*pt*pt + m0*m0*exp(-2.0*lab_rap) + m0*m0*exp(2.0*lab_rap) - 2.0*m0*m0);
		float y_tpc_loedge = log((sqrt(m0*m0+pt*pt*cosh(eta_tpc_loedge)*cosh(eta_tpc_loedge))+pt*sinh(eta_tpc_loedge))/(sqrt(m0*m0+pt*pt)));
		//float y_tpc_hiedge = log((sqrt(m0*m0+pt*pt*cosh(eta_tpc_hiedge)*cosh(eta_tpc_hiedge))+pt*sinh(eta_tpc_hiedge))/(sqrt(m0*m0+pt*pt)));
		float y_tpc_hiedge = 0.0;//same as above
		if( (pt > 0.06) && (lab_rap > y_tpc_loedge) && (lab_rap < y_tpc_hiedge) && (apid == 211 || apid == 321)) {
			refmult3[0]++;
			refmult3[1]++;
		}
                //==================
                if(pid[k] == 2212){//} && y<0. && y>-0.5 ){
                    
                     if(pt > 2.0  || pt < 0.4 ) continue;
                     if(y  > 0.0  || y  < -0.5) continue;
                     //here plot proton accrptance just to test the acceptance is right
                     hpT_y -> Fill( y, pt);
                     np[0]++;
	 	     np[1]++;
		     np[2]++;
		     np[3]++;
                     np[4]++;
                     np[5]++;
                     np[6]++;
                     np[7]++;
                     np[8]++;
                     np[9]++;
                     np[10]++;
                     np[11]++;
                     np[12]++;
                     np[13]++;
                     np[14]++;
                     np[15]++;
                     np[16]++;
                     np[17]++;
                     np[18]++;
                     np[19]++;
                     np[20]++;
                     np[21]++;
                     np[22]++;
                     np[23]++;
                     np[24]++;
                     np[25]++;
                     np[26]++;
		     acceptanceFraction->Fill(40);
                     if(fabs(lab_p)<1.1){np4++; acceptanceFraction->Fill(4);}
                     if(fabs(lab_p)<1.2){np5++; acceptanceFraction->Fill(5);}
                     if(fabs(lab_p)<1.3){np6++; acceptanceFraction->Fill(6);}
		     if(fabs(lab_p)<1.4){np7++; acceptanceFraction->Fill(7);}
                     if(fabs(lab_p)<1.5){np8++; acceptanceFraction->Fill(8);}
                     if(fabs(lab_p)<1.6){np9++; acceptanceFraction->Fill(9);}
                     if(fabs(lab_p)<1.7){np10++; acceptanceFraction->Fill(10);}
                     if(fabs(lab_p)<1.8){np11++; acceptanceFraction->Fill(11);}
                     if(fabs(lab_p)<1.9){np12++;hpT_y_pileuptest -> Fill( y, pt); acceptanceFraction->Fill(12);}
                     if(fabs(lab_p)<2.0){np13++; acceptanceFraction->Fill(13);}
                     if(fabs(lab_p)<2.2){np14++; acceptanceFraction->Fill(14);}
                     if(fabs(lab_p)<2.4){np15++; acceptanceFraction->Fill(15);}
                     if(fabs(lab_p)<2.6){np16++; acceptanceFraction->Fill(16);}
                     if(fabs(lab_p)<2.8){np17++; acceptanceFraction->Fill(17);}
                     if(fabs(lab_p)<3.0){np18++; acceptanceFraction->Fill(18);}
                     if(fabs(lab_p)<3.2){np19++; acceptanceFraction->Fill(19);}
                     if(fabs(lab_p)<3.4){np20++; acceptanceFraction->Fill(20);}
                     if(fabs(lab_p)<3.6){np21++; acceptanceFraction->Fill(21);}
                     if(fabs(lab_p)<3.8){np22++; acceptanceFraction->Fill(22);}
                     if(fabs(lab_p)<4.0){np23++; acceptanceFraction->Fill(23);}
                     if(fabs(lab_p)<4.2){np24++; acceptanceFraction->Fill(24);}
                     if(fabs(lab_p)<4.4){np25++; acceptanceFraction->Fill(25);}
                     if(fabs(lab_p)<4.6){np26++; acceptanceFraction->Fill(26);}
                }
            } //event loop ends
	    if(isPileup){
		refmult3[1]+=pileupRefmult3;
		np[0]+=pileupNp[0];
		np[1]+=pileupNp[1];
		np[2]+=pileupNp[2];
		np[3]+=pileupNp[3];
		np[4]+=pileupNp[4];
                np[5]+=pileupNp[5];
                np[6]+=pileupNp[6];
                np[7]+=pileupNp[7];
                np[8]+=pileupNp[8];
                np[9]+=pileupNp[9];
                np[10]+=pileupNp[10];
                np[11]+=pileupNp[11];
                np[12]+=pileupNp[12];
                np[13]+=pileupNp[13];
                np[14]+=pileupNp[14];
                np[15]+=pileupNp[15];
                np[16]+=pileupNp[16];
                np[17]+=pileupNp[17];
                np[18]+=pileupNp[18];
                np[19]+=pileupNp[19];
                np[20]+=pileupNp[20];
                np[21]+=pileupNp[21];
                np[22]+=pileupNp[22];
                np[23]+=pileupNp[23];
                np[24]+=pileupNp[24];
                np[25]+=pileupNp[25];
                np[26]+=pileupNp[26];
	    }

	    if(refmult3[0]<210 && refmult3[0]>139) npAA->Fill(np[0]);
            if(refmult3[0]<210 && refmult3[0]>139) npAB->Fill(np[1]);
            if(refmult3[1]<210 && refmult3[1]>139) npBA->Fill(np[0]);
            if(refmult3[1]<210 && refmult3[1]>139) npBB->Fill(np[1]);
            otree->Fill();
	    refMult3[0]->Fill(refmult3[0]);
            refMult3[1]->Fill(refmult3[1]);
	    refMult3AndProtons_all->Fill(refmult3[0],np[0]);
            refMult3AndProtons_allsinglesingle->Fill(refmult3[0],np[0]);
            refMult3AndProtons_allsingledouble->Fill(refmult3[0],np[1]);
            refMult3AndProtons_alldoublesingle->Fill(refmult3[1],np[0]);
            refMult3AndProtons_alldoubledouble->Fill(refmult3[1],np[1]);
	    if(isPileup){
	      refMult3AndProtons_singlesingle->Fill(refmult3[0],np[0]);
              refMult3AndProtons_singledouble->Fill(refmult3[0],np[1]);
              refMult3AndProtons_doublesingle->Fill(refmult3[1],np[0]);
              refMult3AndProtons_doubledouble->Fill(refmult3[1],np[1]);
	    }
	    impactParameter->Fill(b);
	    pileupRefmult3=refmult3[0];//If the next event is pileup, just add this value
	    pileupNp[0]=0;
	    pileupNp[1]=np[0];
	    pileupNp[2]=0;
	    pileupNp[3]=np[0];
            pileupNp[4]=np4;
            pileupNp[5]=np5;
            pileupNp[6]=np6;
            pileupNp[7]=np7;
            pileupNp[8]=np8;
            pileupNp[9]=np9;
            pileupNp[10]=np10;
            pileupNp[11]=np11;
            pileupNp[12]=np12;
            pileupNp[13]=np13;
            pileupNp[14]=np14;
            pileupNp[15]=np15;
            pileupNp[16]=np16;
            pileupNp[17]=np17;
            pileupNp[18]=np18;
            pileupNp[19]=np19;
            pileupNp[20]=np20;
            pileupNp[21]=np21;
            pileupNp[22]=np22;
            pileupNp[23]=np23;
            pileupNp[24]=np24;
            pileupNp[25]=np25;
            pileupNp[26]=np26;
        }
        file->Close();
    }// root file loop ends

    output->cd();
    otree->Write();
    TVectorD* v = new TVectorD(1);
    v[0] = pileupFraction;
    v->Write("pileupFraction");
    hpT_y -> Write();
    hpT_y_pileuptest -> Write();
    acceptanceFraction->Write();
    refMult3[0]->Write();
    refMult3[1]->Write();
    refMult3AndProtons_all->Write();
    refMult3AndProtons_allsinglesingle->Write();
    refMult3AndProtons_allsingledouble->Write();
    refMult3AndProtons_alldoublesingle->Write();
    refMult3AndProtons_alldoubledouble->Write();
    refMult3AndProtons_singlesingle->Write();
    refMult3AndProtons_singledouble->Write();
    refMult3AndProtons_doublesingle->Write();
    refMult3AndProtons_doubledouble->Write();
    impactParameter->Write();
    npBB->Write();
    npBA->Write();
    npAB->Write();
    npAA->Write();
    output->Close();
    
    return 0;
}

