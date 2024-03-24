#ifndef LC2
#define LC2

#include "func.h"
#include "algebra.h"

extern int convolution_k;
extern int deduct_signal;
extern int err_bound;
extern int pv_bound;
extern int omit_on;
extern int DOUBLE_PRE;
extern int float_precision;
extern int numeric_interval; 
extern mpf_class empty_mpf;
extern mpf_class* default_mpf;

void Double2MPF (const double&, mpf_class&, int = DOUBLE_PRE);

class LinChiSquare{
	private:
		int dim;
		int simple;
		double max_e;
		double coef_sum;
		int odd;
		int c_mag;
		double o; // copy of odd_c
		mpf_class odd_c; // not used currently
		mpf_class C;
		std::vector< mpf_class > A;
		std::vector< mpf_class > B;
		std::vector< mpf_class > Ci;
		std::vector<int> M;
		std::vector<class PolyUniExp* > BiChi;
		class PolyUniExp* pdf;
		class PolyUniExp* cdf;
		class PolyUniExp* fdf;
	public:
		LinChiSquare(std::vector<double>);
		~LinChiSquare();
		std::string Output(int = 0);
		void Set_M (int);
		void Set_M (std::vector<int>);
		void Use_Default_M(double, int);
		void Compute_Density(void);
		void Single_PDF(int, class PolyDuoExp*);
		void Bi_PDF(class PolyDuoExp*, class PolyDuoExp*);
		double PDF(mpf_class&);
		double CDF(mpf_class&);
		double FDF(mpf_class&);
		void FDF_Exp_Order(std::vector<int>&);
		double Max_Error(void);
		double FDF_Odd(mpf_class&);
		double FDF_Part(double, int, double, double);
		double Compute_Next_Delta(mpf_class& x);
		double Compute_Err_Bound(mpf_class& x);
		int WarnM (void);
};

#endif

