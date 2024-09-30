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
#include "TRandom3.h"
#include "stdlib.h"
#include "EnergyConfig.h"

using namespace std;

int main()
{
    TRandom3 *rand = new TRandom3(0);
    
    int np[2];

    TFile *output = new TFile("outdir/Figure2Sample.root","RECREATE");
    TTree *otree = new TTree("tree","UrQMD Event Tree");
    otree->Branch("np",       np,       "np[2]/I");
    
    TF1* gaus1 = new TF1("gaus1","gaus",-40,120);
    TF1* gaus2 = new TF1("gaus2","gaus",-40,120);
    gaus1->SetParameters(1,40,5);
    gaus2->SetParameters(1,40,15);
    TH1F *gausProtons = new TH1F("gausProtons","gausProtons",100,0,100);
    TH1F *kurtoticProtons = new TH1F("kurtoticProtons","kurtoticProtons",100,0,100);
    
    cout<<"Starting .."<<endl;
    int Entries = 500000;
    for (int j = 0; j< Entries; j++) {
      for(int jm=0;jm<2;++jm){
        np[jm] = 0;
      }
      double single = gaus1->GetRandom();
      np[0]=single;
      np[1]=single;
      if(j%1000==0)np[1]=gaus2->GetRandom();
      gausProtons->Fill(np[0]);
      kurtoticProtons->Fill(np[1]);

      otree->Fill();
    }

    output->cd();
    otree->Write();
    gausProtons->Write();
    kurtoticProtons->Write();
    output->Close();
    
    return 0;
}

