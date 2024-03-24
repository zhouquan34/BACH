#ifndef LC2_ALGEBRA
#define LC2_ALGEBRA

#include "func.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <gmpxx.h>
#include <algorithm>

extern const int default_mag;
extern const int default_omit_mag_bound;
extern const int default_omit_min;
extern const double default_exp_zero;  // remember the input is at most double. so 1e-16 is enough.

extern int default_pv_prec;
extern int omit_x_deg;

class MultiPoly{   //eg.  1+x+0.5x^2+10x^3
	friend class DuoPoly;
	friend class PolyUniExp;
	private: 
		int dmax;  
		std::vector <mpf_class> coef; 
	public:
		MultiPoly(int=0);
		void Copy(const class MultiPoly*);
		void Print(int=0,int=120);
		void Add_Term(int, const mpf_class&);	
		void Add_MP(const class MultiPoly*);
		void Delete_Term(int);
		void Multiply_Scalar(const mpf_class&);
		void Multiply_MP(const class MultiPoly*); 
		void Convert_DP(const class DuoPoly*);
		void Binomial_Expand(class DuoPoly*); // save binomial expansion to a DuoPoly
		void Taylor_I0(int, const mpf_class&); // set the current MultiPoly to be the Taylor expansion of some I0
		int If_Zero(void);
		int Magnitude(void);
		void Compute_Value(const mpf_class&, mpf_class&);
};

class DuoPoly{ //eg. (1+p)+(2+p^2+0.5p^3)x+(3+p+3p^5)x^3
	friend class MultiPoly;
	private:	
		int dmax;
		std::vector< class MultiPoly* > coef;
	public:
		DuoPoly(int=0);
		~DuoPoly();
		void Copy(const class DuoPoly*);
		void Print(int=0,int=120);
		void Print_Dim(void);
		void Add_Term(int, int, const mpf_class&);
		void Add_MP(int, const class MultiPoly*);
		void Add_DP(const class DuoPoly*);
		void Multiply_Scalar(const mpf_class&);
		void Multiply_MP_X(const class MultiPoly*);
		void Multiply_DP(const class DuoPoly*);
		void Convert_MP(const class MultiPoly*);
		void Binomial_Expand(void);
		void Gamma_Int(const mpf_class&, class DuoPoly*, class DuoPoly*);
		void Power_Int(class DuoPoly*);
		void Transform_MP(class MultiPoly*);
		int If_Zero(void);
		void F2X(void);
		void Omit(void);
		int Max_Deg(void);
};

class DuoExpPoly{ //exp(ecoef*x)*(epoly) 
	friend class PolyDuoExp;
	friend class PolyUniExp;
	private: 	
		mpf_class xcoef;
		mpf_class pcoef;
		class DuoPoly* epoly;  
	public:
		DuoExpPoly();
		DuoExpPoly(const mpf_class&);
		DuoExpPoly(const mpf_class&, const mpf_class&);
		~DuoExpPoly();
		void Initialize(void);
		void Copy(const class DuoExpPoly*);
		void Print(int=0,int=120);			
		void Add_DP(const class DuoPoly*);
		void Multiply_Exp_X (const mpf_class&);
		void Multiply_MP_X (const class MultiPoly*);
		void Binomial_Expand(void);
		void Integration(class PolyDuoExp*);
};

class PolyDuoExp{// dep+dep+dep
	friend class PolyUniExp;
	private:
		int term;
		std::vector<class DuoExpPoly*> dexp;
	public:
		PolyDuoExp(); // start with a constant term
		~PolyDuoExp();
		void Copy(const class PolyDuoExp*);		
		void Print(int=0,int=120);
		void Print_Degree(void);
		void Add_Term(const mpf_class&);
		void Add_Term(const mpf_class&, const mpf_class&); // add a term:  exp(-2x)
		void Add_DEP(class DuoExpPoly*);
		void Add_Const_DP(const class DuoPoly*);				
		void Multiply_Exp_X(const mpf_class&);
		void Multiply_MP_X(const class MultiPoly*);
		void Multiply_PDE(const class PolyDuoExp*);
		void Binomial_Expand(void);
		int Combine(void);
		void Integration(void);
		int Max_Deg(void);
		void F2X(void);
		void Omit(void);
};	

class PolyUniExp{
	private:
		int term;
		std::vector<mpf_class> ecoef;
		std::vector<class MultiPoly*> coef;
	public:
		PolyUniExp();
		~PolyUniExp();
		void Print(int=0,int=120);
		void Convert_PDE(class PolyDuoExp*);
		void Integration(class PolyUniExp*);
		void Complementary(class PolyUniExp*, const mpf_class&);
		std::vector<double> Compute_Value(const mpf_class&, mpf_class& v, int, int = 0);
		double Constant(void);
		void Get_Exp_Order(std::vector<mpf_class>& a, std::vector<int>& o);
};




#endif

