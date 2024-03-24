#ifndef LC2_FUNC
#define LC2_FUNC

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <gmpxx.h>
#include <sstream>
#include "gsl_normal.h"

extern const mpf_class PI;
extern const int default_my_exp; // relative error bound = 1e-40, maximum for computing my natural power.
extern const double default_exp_eb;  // default err bound. Used for my own exp
extern const int default_max_precision;  // 1e-80 should be enough for daily purposes, right?

typedef std::vector<int>::size_type    vint_s;
typedef std::vector<double>::size_type    vdb_s;
typedef std::vector<mpf_class>::size_type vmpf_s;
typedef std::vector< std::pair< std::pair<int,int>, double> >::size_type vp1_s;

extern int input_err;
extern int warning;
extern int warn_signal;
extern int verbose;
extern int input_m_max;
extern int fact_max;
extern int input_bessel_n; 
extern std::vector<int> input_ms;
extern std::vector<mpz_class> factorial;
extern std::map<int, std::vector<double> > bessel_taylor_m_cut;
extern std::vector <int> accuracy;
extern std::ofstream LOG;

template <typename T>
void VecPrint (std::vector<T>& v){
	typename std::vector<T>::size_type i; 
	for ( i = 0; i<v.size(); i++){
		LOG << v[i] << " ";	
	}	
	LOG << std::endl;
}

double Max_DB (std::vector<double>&);
void Compute_E(void);
bool CMP_2DB (std::pair< int, std::vector<double> > , std::pair<int, std::vector <double> > );
bool CMP_PairDB(std::pair< std::pair<int,int>, double > , std::pair< std::pair<int,int>, double > );
double ChiSq_One_PV(double);
double ChiSq_One(double);
void Precompute_Fact(void);
int I0(double, int, mpf_class&, double);
void BinCoef(int, int, mpz_class& );
void Compute_Default_M(void);
int Get_Default_M(double, int);
int Compute_I0_Err (mpf_class&, int, mpf_class& ); // return 1 means approx. 0 means exact.
void MinDistance(std::vector<double>& ,  std::vector<std::pair <int, int> >& );
void Compute_Next_Error (mpf_class,  mpf_class, mpf_class, mpf_class, int, mpf_class&);
void Compute_Error_Bound (mpf_class,  mpf_class, mpf_class, mpf_class, int, mpf_class&);
void Exponential (mpf_class, mpf_class&, double = 1e-16); // absolute precision
void DBessel(int);
double Compute_DB (double, int);
void AE(mpf_class&, mpf_class&);

#endif

