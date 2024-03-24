#ifndef LC2_PRE
#define LC2_PRE

#include "lc2.h"

extern const int default_min_err_b; // this is the min relative error bound (10^x) of each term. Simply raise this when you need high precision guaranteed.
extern const int default_pre;
extern const double default_warn_err;

extern int not_compute_p;
extern int better_bound;
extern int user_precision;
extern int input_fact_max;
extern double small_coef;
extern double small_r;
extern std::vector <LinChiSquare*> even_chi;
extern std::string in_file;
extern std::string out_file;
extern std::string log_file;
extern int cout_digit;
extern int next_err;
extern int limit_err;
extern int max_p;
extern double input_x;
extern double x_for_m;
extern int input_omit_deg;
extern std::vector <double> chi_coef;
extern std::vector <double> cmd_xs;
//extern std::ifstream IN;
//extern std::ofstream OUT;

void Upper_P (std::vector<double>&,  class LinChiSquare*&);
void Lower_P (std::vector<double>&,  class LinChiSquare*&);
void Compute_Even_ChiSq(int);
int Scale_Coef(std::vector<double>&, double&, double&);
int P_Bound(std::vector<double>, mpf_class, std::vector<double>&); 
void Better_Bound(std::vector<double>, mpf_class, std::vector<double>&); 
int Estimate_P(double, double); 
void Optimal_M(std::vector<double>&, mpf_class&, const int&, std::vector<int>&); 
void Decide_M (class LinChiSquare*, mpf_class&, std::vector<double>&, std::vector<double> &);
void Start(void);
void End(double, int);
int Read_Para(int, char **);
void List_Para(void);
#endif

