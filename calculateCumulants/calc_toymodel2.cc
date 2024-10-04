#include<iostream>
#include<fstream>
#include<string>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVectorD.h"
#include "TTree.h"
#include "TProfile.h"
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include "Fluc.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "../generateTrees/EnergyConfig.h"

using namespace std;
int main(int argc, char* argv[])
{
    int DATASETINDEX=-1;
    if(ENERGY==3.2)DATASETINDEX=0;
    else if(ENERGY==3.5)DATASETINDEX=1;
    else if(ENERGY==3.9)DATASETINDEX=2;
    else if(ENERGY==4.5)DATASETINDEX=3;
    else if(ENERGY==5.2)DATASETINDEX=4;
    else if(ENERGY==6.2)DATASETINDEX=5;
    else if(ENERGY==7.2)DATASETINDEX=6;
    else if(ENERGY==7.7)DATASETINDEX=7;    

    Int_t INDEX1 = atoi(argv[1]);
    bool DEBUG = true;
    ifstream ins("file.list");
    string line;
    vector<string> vec;
    while(ins>>line){
        vec.push_back(line);
    }
    cout<<"The first file analyzed is "<<vec[0]<<endl;
    
    TFile *inf;
    TTree *tree;
    double acceptanceFrac = 0.0;
    double failureRate = 0.0;
    const int MUL = 350;
    TH1D* DHist[MUL];
    for(int iSafe =0; iSafe<MUL; ++iSafe) {
        DHist[iSafe] = new TH1D(Form("obj%d",iSafe),"", 300, -0.5, 299.5);
    }
    TH1D* DHist_npart[400];
    for(int iSafe =0; iSafe<400; ++iSafe) {
        DHist_npart[iSafe] = new TH1D(Form("obj%d_npart",iSafe),"", 300, -0.5, 299.5);
    }
    TH1D* nProtons[7];
    for(int i=0; i<7; i++) nProtons[i] = new TH1D(Form("nProtons_%d",i),Form("nProtons_%d",i),300,-0.5,299.5);
    TProfile* ANpart = new TProfile("ANpart","",MUL, -0.5, MUL-0.5);
    
    Int_t CentAll[8][7] = {
      { 87, 72,   51,  36,  25,  16, 11},
      {110, 91,   64,  45,  31,  20, 13},
      {139, 115,  81,  56,  38,  25, 15},
      {171, 142,  99,  68,  46,  30, 18},
      {197, 164, 114,  79,  53,  34, 20},
      {0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0}
    };
    Int_t EndBinAll[8] = { 140, 170, 210, 250, 280,0,0,0};    
    Int_t Cent[7] = {CentAll[DATASETINDEX][0],CentAll[DATASETINDEX][1],CentAll[DATASETINDEX][2],CentAll[DATASETINDEX][3],CentAll[DATASETINDEX][4],CentAll[DATASETINDEX][5],CentAll[DATASETINDEX][6]};

    cout<<"Using centrality for "<<ENERGY<<", Data set index: "<<DATASETINDEX<<Form("  : %d, %d, %d, %d, %d, %d, %d",Cent[6],Cent[5],Cent[4],Cent[3],Cent[2],Cent[1],Cent[0])<<endl;
    
    for(unsigned int iLine=0; iLine<vec.size(); ++iLine){
        
        inf = new TFile(vec[iLine].c_str());
        tree = (TTree*)inf->Get("tree");
        Int_t refmult3[2]  = {0};
        Int_t Npart     = 0;
        float b         = 0;
        Int_t np[29] = {0};
	TVectorD *failureRate_pointer = (TVectorD*)inf->Get("failureFraction");
        failureRate = failureRate_pointer->operator()(0);
	TH1F* acceptanceFraction = (TH1F*)inf->Get("acceptanceFraction");
	int thisAcceptanceBin = acceptanceFraction->FindBin(INDEX1);
        int overallAcceptanceBin = acceptanceFraction->FindBin(40);
	acceptanceFrac = (1.0*acceptanceFraction->GetBinContent(thisAcceptanceBin))/(1.0*acceptanceFraction->GetBinContent(overallAcceptanceBin));
	cout<<"acceptance fraction  = "<<acceptanceFrac<<endl;
	        
        tree->SetBranchAddress("refmult3",  &refmult3);
        tree->SetBranchAddress("Npart",       &Npart);
        tree->SetBranchAddress("b",               &b);
        tree->SetBranchAddress("np",              np);
        
        Int_t Entries = tree->GetEntries();
        if(DEBUG) cout<<"Reading Tree .."<<endl;
        for(int nEntry=0; nEntry<Entries; ++nEntry) {
            tree->GetEntry(nEntry);
            if(nEntry !=0 && nEntry%1000000==0) cout<<nEntry<<endl;
	    int thisrefmult3=0;
	    if(INDEX1==0 || INDEX1==1)thisrefmult3=refmult3[0];
            else if(INDEX1>=2)thisrefmult3=refmult3[1];
            int Np = 0;
            Np = np[INDEX1];

            DHist[thisrefmult3]->Fill(Np);
	    DHist_npart[Npart]->Fill(Np);
            ANpart->Fill(thisrefmult3, Npart);
	    for(int i=0; i<7; i++){if(thisrefmult3>Cent[i]){nProtons[i]->Fill(Np);}}
        }
    }
    if(DEBUG) cout<<"Reading Tree Ends"<<endl;
    
    //centrality loop starts
    Int_t StartBin, EndBin = 0;
    const Int_t nObj=17;
    TGraphErrors *OutGraph[nObj];
    const char* cumNames[nObj]= {
        "C1","C2","C3","C4",
        "C5","C6",
        "R21","R32","R42",
        "R51","R62",
        "k1","k2","k3","k4","k5","k6"
    };
    
    for(int i=0; i<nObj; ++i) {
        OutGraph[i] = new TGraphErrors(7);
        OutGraph[i]->SetName(Form("%s",cumNames[i]));
    }
    
    for(int cbin=0; cbin<7; ++cbin) {
        
        if(DEBUG) cout<<"Centrality Bin : "<<cbin<<endl;
        
        StartBin = Cent[cbin];
        if(cbin == 0)
		EndBin=EndBinAll[DATASETINDEX];
        else EndBin = Cent[cbin-1];
        
        //Data point
        Double_t c1,c2,c3,c4,c5,c6,r21,r32,r42,r51,r62;
        c1 = c2 = c3 = c4 = r21 = r32 = r42 = 0;
        Double_t k1,k2,k3,k4,k5,k6;
        k1=k2=k3=k4=k5=k6=0;
        Double_t ek1,ek2,ek3,ek4,ek5,ek6;
        ek1=ek2=ek3=ek4=ek5=ek6=0;
        c5=c6=r51=r62=0;
        Double_t ec1,ec2,ec3,ec4,ec5,ec6=0;
        Double_t er21,er32,er42,er51,er62=0;
        ec1=ec2=ec3=ec4=ec5=ec6=0;
        er21=er32=er42=er51=er62=0;
        Int_t nEvent = 0;
        Int_t AllEvent=0;
        Double_t anpart = 0;
        
        Double_t m[15];
        Double_t u[15];
        for(int iu=0;iu<15;++iu) u[iu]=0;
        for(int im=0;im<15;++im) m[im]=0;
        
        if(DEBUG) cout<<"Data point loop start .."<<endl;
        for(int ibin=StartBin; ibin<EndBin; ++ibin) {
            
            //nEvent = DHist_npart[ibin]->Integral();
	    nEvent = DHist[ibin]->Integral();
            if(nEvent<1) continue;
            AllEvent += nEvent;
            Fluc fluc;
            fluc.SetMaxOrder(12);
            fluc.Init();
            //fluc.ReadHistogram(DHist_npart[ibin]);
            fluc.ReadHistogram(DHist[ibin]);
            for(int im=1;im<13;++im){
                m[im] = fluc.GetMoment(im);
            }
            if(m[1] == 0 ) continue;
            if(isnan(m[1])) continue;
            for(int iu=2;iu<13;++iu){
                u[iu] = fluc.GetCentralMoment(iu);
            }
            Double_t sigma = sqrt(u[2]);
            if(isnan(sigma)) {
                continue;
            }
            if(sigma == 0) continue;
            
            Double_t dc1 = fluc.GetCumulant(1);
            Double_t dc2 = fluc.GetCumulant(2);
            Double_t dc3 = fluc.GetCumulant(3);
            Double_t dc4 = fluc.GetCumulant(4);
            Double_t dc5 = fluc.GetCumulant(5);
            Double_t dc6 = fluc.GetCumulant(6);
            
            Double_t dk2 = dc2 - dc1;
            Double_t dk3 = 2*dc1 - 3*dc2 + dc3;
            Double_t dk4 = -6*dc1 + 11*dc2 - 6*dc3 + dc4;
            Double_t dk5 = 24*dc1 - 50*dc2 + 35*dc3 - 10*dc4 + dc5;
            Double_t dk6 = -120*dc1 + 274*dc2 - 225*dc3 + 85*dc4 - 15*dc5 + dc6;
            
            k2 += nEvent * dk2;
            k3 += nEvent * dk3;
            k4 += nEvent * dk4;
            k5 += nEvent * dk5;
            k6 += nEvent * dk6;
            
            c1  += nEvent * dc1;
            c2  += nEvent * dc2;
            c3  += nEvent * dc3;
            c4  += nEvent * dc4;
            c5  += nEvent * dc5;
            c6  += nEvent * dc6;
            
            if(dc1!=0) r21 += nEvent * (dc2/dc1);
            if(dc2!=0) r32 += nEvent * (dc3/dc2);
            if(dc2!=0) r42 += nEvent * (dc4/dc2);
            if(dc1!=0) r51 += nEvent * (dc5/dc1);
            if(dc2!=0) r62 += nEvent * (dc6/dc2);
            
            ec1 += nEvent * fabs(u[2]);
            ec2 += nEvent * fabs(u[4]-u[2]*u[2]);
            ec3 += nEvent * fabs(u[6]-u[3]*u[3]-6*u[4]*u[2]+9*pow(u[2],3));
            ec4 += nEvent * fabs((-36*pow(u[2], 4) + 48*pow(u[2], 2)*u[4] + 64*u[2]*pow(u[3], 2) - 12*u[2]*u[6] - 8*u[3]*u[5] - pow(u[4], 2) + u[8]));
            ec5 += nEvent * fabs(u[10] + 900*pow(u[2],5) - 900*pow(u[2],3)*u[4] - 1000*pow(u[2],2)*pow(u[3],2) + 160*pow(u[2],2)*u[6] + 240*u[2]*u[3]*u[5] + 125*u[2]*pow(u[4],2) - 20*u[2]*u[8] + 200*pow(u[3],2)*u[4] - 20*u[3]*u[7] - 10*u[4]*u[6] + pow(u[5],2));
            ec6 += nEvent * fabs(-30*u[10]*u[2] + u[12] - 8100*pow(u[2],6) + 13500*pow(u[2],4)*u[4] + 39600*pow(u[2],3)*pow(u[3],2) - 2880*pow(u[2],3)*u[6] - 9720*pow(u[2],2)*u[3]*u[5] - 3600*pow(u[2],2)*pow(u[4],2) + 405*pow(u[2],2)*u[8] - 9600*u[2]*pow(u[3],2)*u[4] + 840*u[2]*u[3]*u[7] + 510*u[2]*u[4]*u[6] + 216*u[2]*pow(u[5],2) - 400*pow(u[3],4) + 440*pow(u[3],2)*u[6] + 1020*u[3]*u[4]*u[5] - 40*u[3]*u[9] + 225*pow(u[4],3) - 30*u[4]*u[8] - 12*u[5]*u[7] - pow(u[6],2));
            
            er21 += nEvent * fabs(-pow(m[1],2) + m[2] - 2*m[3]/m[1] + 2*pow(m[2],2)/pow(m[1],2) + m[4]/pow(m[1],2) - 2*m[2]*m[3]/pow(m[1],3) + pow(m[2],3)/pow(m[1],4));
            
            er32 += nEvent * fabs( (9*u[2] - 6*u[4]/u[2] + 6*pow(u[3], 2)/pow(u[2], 2) + u[6]/pow(u[2], 2) - 2*u[3]*u[5]/pow(u[2], 3) + pow(u[3], 2)*u[4]/pow(u[2], 4)) );
            
            er42 += nEvent * fabs((-9*pow(u[2], 2) + 9*u[4] + 40*pow(u[3], 2)/u[2] - 6*u[6]/u[2] - 8*u[3]*u[5]/pow(u[2], 2) + 6*pow(u[4], 2)/pow(u[2], 2) + u[8]/pow(u[2], 2) + 8*pow(u[3], 2)*u[4]/pow(u[2], 3) - 2*u[4]*u[6]/pow(u[2], 3) + pow(u[4], 3)/pow(u[2], 4)));
            
            er51 += nEvent * fabs( (m[10]/pow(m[1], 2) - 9216*pow(m[1], 8) + 43776*pow(m[1], 6)*m[2] - 19200*pow(m[1], 5)*m[3] - 66960*pow(m[1], 4)*pow(m[2], 2) + 8400*pow(m[1], 4)*m[4] + 46080*pow(m[1], 3)*m[2]*m[3] - 3360*pow(m[1], 3)*m[5] + 36000*pow(m[1], 2)*pow(m[2], 3) - 15720*pow(m[1], 2)*m[2]*m[4] - 5920*pow(m[1], 2)*pow(m[3], 2) + 1192*pow(m[1], 2)*m[6] - 25680*m[1]*pow(m[2], 2)*m[3] + 4608*m[1]*m[2]*m[5] + 2400*m[1]*m[3]*m[4] - 320*m[1]*m[7] - 3600*pow(m[2], 4) + 6600*pow(m[2], 2)*m[4] + 4800*m[2]*pow(m[3], 2) - 1240*m[2]*m[6] - 480*m[3]*m[5] - 25*pow(m[4], 2) + 65*m[8] - 1200*pow(m[2], 3)*m[3]/m[1] - 960*pow(m[2], 2)*m[5]/m[1] - 1300*m[2]*m[3]*m[4]/m[1] + 220*m[2]*m[7]/m[1] - 400*pow(m[3], 3)/m[1] + 140*m[3]*m[6]/m[1] - 40*m[4]*m[5]/m[1] - 10*m[9]/m[1] + 1500*pow(m[2], 2)*pow(m[3], 2)/pow(m[1], 2) + 100*pow(m[2], 2)*m[6]/pow(m[1], 2) - 60*m[2]*m[3]*m[5]/pow(m[1], 2) - 20*m[2]*m[8]/pow(m[1], 2) + 100*pow(m[3], 2)*m[4]/pow(m[1], 2) - 20*m[3]*m[7]/pow(m[1], 2) + 10*pow(m[5], 2)/pow(m[1], 2) - 200*pow(m[2], 2)*m[3]*m[4]/pow(m[1], 3) - 200*m[2]*pow(m[3], 3)/pow(m[1], 3) + 20*m[2]*m[3]*m[6]/pow(m[1], 3) + 20*m[2]*m[4]*m[5]/pow(m[1], 3) + 20*pow(m[3], 2)*m[5]/pow(m[1], 3) - 2*m[5]*m[6]/pow(m[1], 3) + 100*pow(m[2], 3)*pow(m[3], 2)/pow(m[1], 4) - 20*pow(m[2], 2)*m[3]*m[5]/pow(m[1], 4) + m[2]*pow(m[5], 2)/pow(m[1], 4))  );
            
            er62 += nEvent * fabs(-30*u[10]/u[2] + u[12]/pow(u[2],2) - 3600*pow(u[2],4) + 5400*pow(u[2],2)*u[4] + 30000*u[2]*pow(u[3],2) - 1800*u[2]*u[6] - 8160*u[3]*u[5] - 225*pow(u[4],2)+ 345*u[8] - 3900*pow(u[3],2)*u[4]/u[2] + 840*u[3]*u[7]/u[2] - 120*u[4]*u[6]/u[2] + 216*pow(u[5],2)/u[2] + 2300*pow(u[3],4)/pow(u[2],2) - 140*pow(u[3],2)*u[6]/pow(u[2],2) + 240*u[3]*u[4]*u[5]/pow(u[2],2) - 40*u[3]*u[9]/pow(u[2],2) - 12*u[5]*u[7]/pow(u[2],2) + 30*pow(u[6],2)/pow(u[2],2) - 520*pow(u[3],3)*u[5]/pow(u[2],3) + 20*pow(u[3],2)*u[8]/pow(u[2],3) + 52*u[3]*u[5]*u[6]/pow(u[2],3) - 2*u[6]*u[8]/pow(u[2],3) + 100*pow(u[3],4)*u[4]/pow(u[2],4) - 20*pow(u[3],2)*u[4]*u[6]/pow(u[2],4) + u[4]*pow(u[6],2)/pow(u[2],4));
            
            ek2 += nEvent * fabs(-pow(u[2], 2) + u[2] - 2*u[3] + u[4]);
            ek3 += nEvent * fabs(9*pow(u[2], 3) - 21*pow(u[2], 2) + 24*u[2]*u[3] - 6*u[2]*u[4] + 4*u[2] - pow(u[3], 2) - 12*u[3] + 13*u[4] - 6*u[5] + u[6]);
            ek4 += nEvent * fabs(-36*pow(u[2], 4) + 456*pow(u[2], 3) - 432*pow(u[2], 2)*u[3] + 48*pow(u[2], 2)*u[4] - 337*pow(u[2], 2) + 64*u[2]*pow(u[3], 2) + 648*u[2]*u[3] - 370*u[2]*u[4] + 108*u[2]*u[5] - 12*u[2]*u[6] + 36*u[2] - 124*pow(u[3], 2) + 60*u[3]*u[4] - 8*u[3]*u[5] - 132*u[3] - pow(u[4], 2) + 193*u[4] - 144*u[5] + 58*u[6] - 12*u[7] + u[8]);
            ek5 += nEvent * fabs(u[10] + 900*pow(u[2], 5) - 9900*pow(u[2], 4) + 8400*pow(u[2], 3)*u[3] - 900*pow(u[2], 3)*u[4] + 18465*pow(u[2], 3) - 1000*pow(u[2], 2)*pow(u[3], 2) - 30200*pow(u[2], 2)*u[3] + 10050*pow(u[2], 2)*u[4] - 1920*pow(u[2], 2)*u[5] + 160*pow(u[2], 2)*u[6] - 7540*pow(u[2], 2) + 9900*u[2]*pow(u[3], 2) - 3400*u[2]*u[3]*u[4] + 240*u[2]*u[3]*u[5] + 18800*u[2]*u[3] + 125*u[2]*pow(u[4], 2) - 15070*u[2]*u[4] + 7400*u[2]*u[5] - 2110*u[2]*u[6] + 320*u[2]*u[7] - 20*u[2]*u[8] + 576*u[2] - 800*pow(u[3], 3) + 200*pow(u[3], 2)*u[4] - 5705*pow(u[3], 2) + 5000*u[3]*u[4] - 1570*u[3]*u[5] + 280*u[3]*u[6] - 20*u[3]*u[7] - 2400*u[3] - 450*pow(u[4], 2) + 120*u[4]*u[5] - 10*u[4]*u[6] + 4180*u[4] - pow(u[5], 2) - 3980*u[5] + 2273*u[6] - 800*u[7] + 170*u[8] - 20*u[9]);
            ek6 += nEvent * fabs(-30*u[10]*u[2] + 395*u[10] - 30*u[11] + u[12] - 8100*pow(u[2], 6) + 294300*pow(u[2], 5) - 243000*pow(u[2], 4)*u[3] + 13500*pow(u[2], 4)*u[4] - 916920*pow(u[2], 4) + 39600*pow(u[2], 3)*pow(u[3], 2) + 1395000*pow(u[2], 3)*u[3] - 340200*pow(u[2], 3)*u[4] + 48600*pow(u[2], 3)*u[5] - 2880*pow(u[2], 3)*u[6] + 843105*pow(u[2], 3) - 510600*pow(u[2], 2)*pow(u[3], 2) + 144000*pow(u[2], 2)*u[3]*u[4] - 9720*pow(u[2], 2)*u[3]*u[5] - 1838400*pow(u[2], 2)*u[3] - 3600*pow(u[2], 2)*pow(u[4], 2) + 918810*pow(u[2], 2)*u[4] - 313650*pow(u[2], 2)*u[5] + 67620*pow(u[2], 2)*u[6] - 8100*pow(u[2], 2)*u[7] + 405*pow(u[2], 2)*u[8] - 237076*pow(u[2], 2) + 48000*u[2]*pow(u[3], 3) - 9600*u[2]*pow(u[3], 2)*u[4] + 876620*u[2]*pow(u[3], 2) - 548250*u[2]*u[3]*u[4] + 115200*u[2]*u[3]*u[5] - 14700*u[2]*u[3]*u[6] + 840*u[2]*u[3]*u[7] + 697200*u[2]*u[3] + 48525*u[2]*pow(u[4], 2) - 10350*u[2]*u[4]*u[5] + 510*u[2]*u[4]*u[6] - 683810*u[2]*u[4] + 216*u[2]*pow(u[5], 2) + 439710*u[2]*u[5] - 183218*u[2]*u[6] + 48900*u[2]*u[7] - 8070*u[2]*u[8] + 750*u[2]*u[9] + 14400*u[2] - 400*pow(u[3], 4) - 111000*pow(u[3], 3) + 72200*pow(u[3], 2)*u[4] - 8400*pow(u[3], 2)*u[5] + 440*pow(u[3], 2)*u[6] - 272945*pow(u[3], 2) - 9750*u[3]*pow(u[4], 2) + 1020*u[3]*u[4]*u[5] + 322950*u[3]*u[4] - 146298*u[3]*u[5] + 45150*u[3]*u[6] - 8580*u[3]*u[7] + 900*u[3]*u[8] - 40*u[3]*u[9] - 65760*u[3] + 225*pow(u[4], 3) - 49195*pow(u[4], 2) + 24750*u[4]*u[5] - 4970*u[4]*u[6] + 600*u[4]*u[7] - 30*u[4]*u[8] + 129076*u[4] - 1245*pow(u[5], 2) + 210*u[5]*u[6] - 12*u[5]*u[7] - 143700*u[5] - pow(u[6], 2) + 100805*u[6] - 46710*u[7] + 14523*u[8] - 3000*u[9]);
            
            //if(nEvent > 0) anpart +=  (ANpart->GetBinContent(ibin+1) * nEvent);
	    if(nEvent > 0) anpart +=  (ANpart->GetBinContent(ibin+1) * nEvent);
            //if(nEvent > 0) anpart +=  ((ibin+1) * nEvent);
        }
        
        c1  /= AllEvent;
        c2  /= AllEvent;
        c3  /= AllEvent;
        c4  /= AllEvent;
        c5  /= AllEvent;
        c6  /= AllEvent;
        
        k1   = c1;
        k2  /= AllEvent;
        k3  /= AllEvent;
        k4  /= AllEvent;
        k5  /= AllEvent;
        k6  /= AllEvent;
        
        r21 = r21/AllEvent;
        r32 = r32/AllEvent;
        r42 = r42/AllEvent;
        r51 = r51/AllEvent;
        r62 = r62/AllEvent;
        
        ec1 = sqrt(ec1) / AllEvent;
        ec2 = sqrt(ec2) / AllEvent;
        ec3 = sqrt(ec3) / AllEvent;
        ec4 = sqrt(ec4) / AllEvent;
        ec5 = sqrt(ec5) / AllEvent;
	ec6 = sqrt(ec6) / AllEvent;
        //c6 = sqrt(ec6) / AllEvent;
        
        ek1 = ec1;
        ek2 = sqrt(ek2) / AllEvent;
        ek3 = sqrt(ek3) / AllEvent;
        ek4 = sqrt(ek4) / AllEvent;
        ek5 = sqrt(ek5) / AllEvent;
        ek6 = sqrt(ek6) / AllEvent;
        
        er21 = sqrt(er21) / AllEvent;
        er32 = sqrt(er32) / AllEvent;
        er42 = sqrt(er42) / AllEvent;
        er51 = sqrt(er51) / AllEvent;
        er62 = sqrt(er62) / AllEvent;
        
        anpart /= AllEvent;
        cout<<"average npart : "<<anpart<<endl;
        
        OutGraph[0]  ->SetPoint(cbin, anpart, c1);
        OutGraph[1]  ->SetPoint(cbin, anpart, c2);
        OutGraph[2]  ->SetPoint(cbin, anpart, c3);
        OutGraph[3]  ->SetPoint(cbin, anpart, c4);
        OutGraph[4]  ->SetPoint(cbin, anpart, c5);
        OutGraph[5]  ->SetPoint(cbin, anpart, c6);
        OutGraph[6]  ->SetPoint(cbin, anpart, r21);
        OutGraph[7]  ->SetPoint(cbin, anpart, r32);
        OutGraph[8]  ->SetPoint(cbin, anpart, r42);
        OutGraph[9]  ->SetPoint(cbin, anpart, r51);
        OutGraph[10] ->SetPoint(cbin, anpart, r62);
        OutGraph[11] ->SetPoint(cbin, anpart, k1);
        OutGraph[12] ->SetPoint(cbin, anpart, k2);
        OutGraph[13] ->SetPoint(cbin, anpart, k3);
        OutGraph[14] ->SetPoint(cbin, anpart, k4);
        OutGraph[15] ->SetPoint(cbin, anpart, k5);
        OutGraph[16] ->SetPoint(cbin, anpart, k6);
        
        OutGraph[0]  ->SetPointError(cbin, 0, ec1);
        OutGraph[1]  ->SetPointError(cbin, 0, ec2);
        OutGraph[2]  ->SetPointError(cbin, 0, ec3);
        OutGraph[3]  ->SetPointError(cbin, 0, ec4);
        OutGraph[4]  ->SetPointError(cbin, 0, ec5);
        OutGraph[5]  ->SetPointError(cbin, 0, ec6);
        OutGraph[6]  ->SetPointError(cbin, 0, er21);
        OutGraph[7]  ->SetPointError(cbin, 0, er32);
        OutGraph[8]  ->SetPointError(cbin, 0, er42);
        OutGraph[9]  ->SetPointError(cbin, 0, er51);
        OutGraph[10] ->SetPointError(cbin, 0, er62);
        OutGraph[11] ->SetPointError(cbin, 0, ek1);
        OutGraph[12] ->SetPointError(cbin, 0, ek2);
        OutGraph[13] ->SetPointError(cbin, 0, ek3);
        OutGraph[14] ->SetPointError(cbin, 0, ek4);
        OutGraph[15] ->SetPointError(cbin, 0, ek5);
        OutGraph[16] ->SetPointError(cbin, 0, ek6);
        //Data point loop ends
        if(DEBUG) cout<<"Data point loop ends .."<<endl;
        
        for(int ifree=StartBin; ifree<EndBin; ++ifree) {
            //delete DHist_npart[ifree];
	    delete DHist[ifree];
        }//delete Histograms of p pbar
        
    }// centrality loop ends
    
    if(DEBUG) cout<<"Writing .."<<endl;
    TString outName = "dummy.root";
    outName = Form("Output/cum_%d_%ffailureRate.root", INDEX1,failureRate);
    
    TFile *out = new TFile(outName.Data(), "recreate");
    
    for(int i=0; i<7; i++) nProtons[i]->Write();
    out->cd();
    for(int i=0; i<17; ++i) {
        OutGraph[i]->Write();
    }
    TVectorD* v = new TVectorD(1);
    v[0] = failureRate;
    v->Write("failureRate");
    TVectorD* v2 = new TVectorD(1);
    v2[0] = acceptanceFrac;
    v2->Write("acceptanceFrac");
    out->Close();
    
    delete out;
    
    if(DEBUG) cout<<"Ending .."<<endl;
    
    return 0;
}

