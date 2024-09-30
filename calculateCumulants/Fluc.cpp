
#include <iostream>
#include "Fluc.h" 

using std::cout;
using std::cerr;
using std::endl;
using std::vector;

//__________________________________________________________________________________
Fluc::Fluc(){
	maxOrder = 4;
}

void Fluc::Reset(){
    maxOrder = 4;
    ncm.clear();
    ncm_cor.clear();
    fm.clear();
    C.clear();
    C_cor.clear();

}
//__________________________________________________________________________________
void Fluc::SetMaxOrder( const int order ){
	maxOrder = order;
}

//__________________________________________________________________________________
//  In order to calculate the error of m-th cumulant, 
//  up to 2m-th cumulants have to be calculated.
//__________________________________________________________________________________
void Fluc::SetDeltaTheoremOn(){
	maxOrder = 2*maxOrder;
}

//__________________________________________________________________________________
void Fluc::Init(){
	ncm.resize(2*maxOrder+1);
	ncm_cor.resize(2*maxOrder+1);
	fm.resize(maxOrder+1);
	C.resize(2*maxOrder+1);
	C_cor.resize(2*maxOrder+1);
}

//__________________________________________________________________________________
//  Read 1-D histogram and non-central moments up to 2*maxorder is calculated. 
//  This function cannot be used for efficiency correction.
//  2-D information is needed for efficiency correction.	
//  Factorial moments are of "volume fluctuation" components.
//__________________________________________________________________________________
void Fluc::ReadHistogram( TH1* h ){
	nevent = h->Integral();
	const int nxbin = h->GetXaxis()->GetNbins();
	for(int ixbin=0; ixbin<nxbin; ixbin++){
		const LongDouble_t nBinEvent = h->GetBinContent(ixbin+1);
		if(nBinEvent==0) continue;
	//	const LongDouble_t m = h->GetBinLowEdge(ixbin+1);
		const LongDouble_t m = h->GetBinCenter(ixbin+1);
		for(int r=0; r<=maxOrder; r++){
			ncm[r] += pow(m,r)*nBinEvent;
		}
	}
	for(int r=0; r<=maxOrder; r++){
		ncm[r] /= nevent;
	}
}


//__________________________________________________________________________________
void Fluc::SetMoments( const vector<LongDouble_t> &vtmpNCM ){
	for(int i=1; i<=maxOrder; i++){
		ncm[i] = vtmpNCM[i];
	}
}

//__________________________________________________________________________________
void Fluc::SetMoments( const int n, const LongDouble_t val ){
	ncm[n] = val;
}

//__________________________________________________________________________________
void Fluc::SetCumulants( const int n, const LongDouble_t val ){
	C[n] = val;
}


//__________________________________________________________________________________
//      Efficiency correction on <m1m2> terms with Stirling numbers
//__________________________________________________________________________________
LongDouble_t Fluc::GetEffm1m2( int a, int b ){
	if(a==0&&b==0) return 1;
	LongDouble_t rterm = 0;
	for(int i=0; i<=a; i++){
	for(int j=0; j<=b; j++){
		if(i==a&&j==b) continue;
		rterm += GetStirlingNumber(a,i)*GetStirlingNumber(b,j)
			*( m1m2[i][j]/(pow(eff1,a)*pow(eff2,b)) - GetEffm1m2(i,j) );
	}
	}
	return m1m2[a][b]/(pow(eff1,a)*pow(eff2,b)) + rterm;
}

//__________________________________________________________________________________
LongDouble_t Fluc::GetCumulant( int n ){
	if(n>maxOrder){
		cerr << n <<"-th cumulant has not been calculated."<< endl;
		return -9999;
	}
	if(n==1) return ncm[1]; 
	C[n] = ncm[n];
	for(int m=1; m<=n-1; m++){
		C[n] -= TMath::Binomial(n-1,m-1)*GetCumulant(m)*ncm[n-m];
	}
	return C[n];
}

//_________________________________________________________________________________________
LongDouble_t Fluc::GetMoment( const int n ) const {
	return ncm[n];
}

//__________________________________________________________________________________
LongDouble_t Fluc::GetMomentsFromCumulants( int n ){
	if(n==0) return 1;
	if(n==1) return C[1];
	LongDouble_t val = 0;
	for(int i=0; i<=n-1; i++){
		val += TMath::Binomial(n-1,i)*C[i+1]*GetMomentsFromCumulants(n-1-i);
	}
	return val;
}

//_________________________________________________________________________________________
LongDouble_t Fluc::GetCentralMoment( const int n ) {
	if(n>maxOrder){
		cerr << n <<"-th moment has not been calculated."<< endl;
		return -9999;
	}
	LongDouble_t cm = 0;
	for(int k=0; k<=n; k++){
		cm += TMath::Binomial(n,k)*pow(-1,n-k)*ncm[k]*pow(ncm[1],n-k);
	}
	return cm;
}

//__________________________________________________________________________________
LongDouble_t Fluc::GetEffCumulant( int n ){
	if(n>maxOrder){
		cerr << n <<"-th cumulant has not been calculated."<< endl;
		return -9999;
	}
	if(n==1) return ncm_cor[1]; 
	C_cor[n] = ncm_cor[n];
	for(int m=1; m<=n-1; m++){
		C_cor[n] -= TMath::Binomial(n-1,m-1)*GetEffCumulant(m)*ncm_cor[n-m];
	}
	return C_cor[n];
}

//_________________________________________________________________________________________
LongDouble_t Fluc::GetFactorialMoment( const int n ){
	LongDouble_t fact = 0;
	for(int i=0; i<=n; i++){
		fact += ncm[i]*(LongDouble_t)GetStirlingNumber(n,i);
	}
	return fact;
}

//__________________________________________________________________________________
int Fluc::StirlingNumber( int n, int k ){
	if(n==k) return 1;
	else if(n==0||k==0) return 0;
	else{
		return StirlingNumber(n-1,k-1) 
			+ (n-1)*StirlingNumber(n-1,k);
	}
}


//__________________________________________________________________________________
int Fluc::GetStirlingNumber( const int n, const int k ){
	if((n-k)%2!=0) return -StirlingNumber(n,k);
	else return StirlingNumber(n,k);
}

