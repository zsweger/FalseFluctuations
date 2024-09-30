#ifndef __Fluc_h__
#define __Fluc_h__

#include <vector>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>

//__________________________________________________________________________________
class Fluc {

	private :
		LongDouble_t nevent;
		LongDouble_t eff1, eff2;
		std::vector<LongDouble_t> ncm;
		std::vector<LongDouble_t> C;
		std::vector< std::vector<LongDouble_t> > m1m2;
		std::vector<LongDouble_t> fm;  // factorial moment for 1D histogram
		std::vector<LongDouble_t> ncm_cor;
		std::vector< std::vector<LongDouble_t> > m1m2_cor;
		std::vector<LongDouble_t> C_cor;

	protected :
		int maxOrder;

	public :
		Fluc();
        void Reset();
        //void ReadList(TFile* f);
		virtual ~Fluc(){};
		void SetMaxOrder( const int order);
		void SetDeltaTheoremOn(); 
		void Init();
		void ReadHistogram( TH1* h );
		// Setter
		void SetMoments( const std::vector<LongDouble_t> &vtmpNCM );
		void SetMoments( const int n, const LongDouble_t val );
		void SetCumulants( const int n, const LongDouble_t val );
		// Getter
		LongDouble_t GetEffm1m2( int a, int b );
		LongDouble_t GetCumulant( int n );
		LongDouble_t GetEffCumulant( int n );
		LongDouble_t GetMoment( const int n ) const;
		LongDouble_t GetCentralMoment( const int n );
		LongDouble_t GetMomentsFromCumulants( int n );
		LongDouble_t GetFactorialMoment( const int n );
		int StirlingNumber( int n, int k );
		int GetStirlingNumber( const int n, const int k );

};

#endif
