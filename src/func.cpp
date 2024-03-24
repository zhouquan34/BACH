/*
   Copyright 09/30/2015 Quan Zhou, Yongtao Guan, Baylor College of Medicine
   This file is part of program bach. You can redistribute it or modify it under the terms of the GNU General Public License. See it for details. You should have received a copy of GPL which is named 'COPYING'.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "func.h"

using namespace std;

const mpf_class PI ("3.14159265358979323846264338327950288419716939937510");
const int default_my_exp = 185;
const double default_exp_eb = 1e-40;
const int default_max_precision = 80;

vector<mpf_class> E;
const int default_e_size = 10;
const int default_e1_degree = 80; // error < 1e-118

ofstream LOG;

vector <mpz_class> factorial;
map < int, vector <double> > bessel_taylor_m_cut;
vector <int> accuracy;
vector <int> input_ms;
int warn_signal = 0;
int warning = 1;
int verbose = 1;
int fact_max = 500;
int input_m_max = 80;
int input_bessel_n = -1;
int input_err = 0;

double Max_DB (vector<double>& x){
	if (x.size() == 0){
		cerr << "Wrong! Empty vector for Max_DB!" << endl;
		return 0.0;
	}
	double max = x[0];
	for (int i=0; i< (int) x.size(); i++){
		if (x[i] > max){max = x[i];}	
	}
	return max;
}

void Compute_E (void){ 	
	mpf_class tmp ("0");
	E.push_back(tmp);
	for (int i=0; i<=default_e1_degree; i++){
		mpf_class t0 (factorial[i]);
		tmp += 1/t0;
	}
	E.push_back(tmp);
	for (int k = 2; k<=default_e_size; k++){
		tmp = tmp*tmp;
		E.push_back(tmp);
	}
}

bool CMP_2DB (pair< int, vector<double> > x, pair<int, vector <double> > y){
	if  (x.second[0] == y.second[0]){
		return x.second[1]>y.second[1];
	}else{
		return x.second[0]>y.second[0];	
	}
}

bool CMP_PairDB(pair< pair<int,int>, double > x, pair< pair<int,int>, double > y){return x.second<y.second;}

double ChiSq_One(double x){
	return (exp(-x/2.0)/sqrt(2.0* (PI.get_d()) *x));	
}

double ChiSq_One_PV(double x){return 2.0*gsl_cdf_ugaussian_P(-sqrt(x));}

void Precompute_Fact (void){
	if (verbose >= 6) LOG << "Computing factorial" <<  endl;
	mpz_class f0(1);
	factorial.push_back(f0);
	for (int i=1; i<=fact_max; i++){
		mpz_class f(factorial.at(i-1));	
		f=f*i;
		factorial.push_back(f);
	}
	return;
} 

int I0 (double x, int n, mpf_class& bi, double err){ // This function is only used in this cpp. So x is double
	if (n > fact_max){n = fact_max;}
	mpf_set_ui(bi.get_mpf_t(), 1);
	int m = 1;
	mpf_class s1 (bi);
	mpf_class hx (x/2);
	mpf_class s (0);
	mpf_class q (1);
	mpf_class q1(1);
//	mpf_class q2 (1);
	mpf_class rmd (0);
	for (m=1; m<=n; m++){
		q = hx * hx / (m*m);
		s = s1 * q;
		if (verbose >= 10) LOG << "m = " << m << "; q = " << q << "; s = " << s << endl;
		if ( q < 1){
//			q1 = hx*hx/( (m+1)*(m+1) );
//			q2 = hx*hx/( (m+2)*(m+2) );
//			rmd = s + s*q1 + s*q1*q2/(1-q2);  // this is a bound 
			q1 = hx * hx /( (m+1)*(m+1) );
			rmd = s/(1-q1);
			if (rmd / bi < err ){
				if (verbose >= 9) LOG << "I0 of " << x << " = " << bi << "; stop at " << m  << "; rmd = " << rmd  << "; q = "<< q<< endl;
				break;
			}
		}
		bi += s;
		s1 = s;
	}
	return (m-1);
}

int Compute_I0_Err (mpf_class& x, int n, mpf_class& e){
	int m = 1;
	mpf_class bi (1);
	mpf_class s1 (bi);
	mpf_class hx (x/2);
	mpf_class s (0);
	mpf_class q (1);
	mpf_class q1 (1);
	mpf_class rmd (0);
	for (m=1; m<=n; m++){
		q = hx * hx / (m*m);
		s = s1 * q;
		bi += s;
		s1 = s;
	}
	q = hx*hx/((n+1)*(n+1));
	s = s1*q;
	q1 = hx*hx/((n+2)*(n+2));
	if (verbose >=4) LOG << "In Compute_I0_Err: x = " << x << "; n = " << n << "; bi = "  << bi << endl;
	if (q1 < 1){
		rmd = s/(1-q1);
		e = rmd/bi;
		if (verbose >=4) LOG <<  "rmd = " << rmd << "; e =" <<  e<< "; q=" << q<<  endl;
		return 0;
	}else{
		s1 = s;
		rmd = s;
		for (m=n+2; m<=n+10; m++){
			q = hx * hx / (m*m);
			s = s1 * q;
			rmd += s;
			s1 = s;
		}	
		e = rmd/bi;
		if (verbose >=4) LOG <<  "rmd = " << rmd << "; e ="<< e << endl;
		return 1;
	}
}


void BinCoef(int n, int k, mpz_class& coef){
	coef = factorial[n]/( factorial[k] * factorial[n-k]);
	if (verbose >= 10) LOG << "Computing coef of " << n << "C" << k  << " = " << coef << endl;
	return;
}

void Compute_Default_M (void){
	if (input_err == 0){
		for (int i=1; i<=10; i++) accuracy.push_back(i);
		for (int i=12; i<=default_max_precision; i+=4) accuracy.push_back(i);
	}else{
		accuracy.push_back(input_err);		
	}
//	if (verbose >= 7) LOG << "acc size = " << accuracy.size() << endl;

	for (int i=0; i< (int) accuracy.size(); i++){
		if (verbose >= 7) LOG << "At acc = " << accuracy[i] << "; i = " << i << endl;
		vector <double> bt_cut (input_m_max+1, 0.0);
		double eb = pow(10.0, -accuracy[i]);
		mpf_class tmp;
		for (int z= -10; z<= -1; z++){
			double x = pow(10.0, z);
			int m = I0(x, fact_max, tmp, eb);
			if (m > input_m_max) {bt_cut[input_m_max] = x;}
			else { bt_cut[m] = x; }
		}
		int z = 2;
		while (1){
//			LOG << "Z = " << z << endl;
			double x = ((double) z) * 0.1;
			int m = I0(x, fact_max, tmp, eb);
			if (m>input_m_max){break;}
			bt_cut[m] = x;
			z ++;
		}
		bessel_taylor_m_cut[accuracy[i]] = bt_cut;
		if (verbose >= 7) {
			LOG << "Default M selection setup for error bound = "<< eb  << endl;
//			LOG << input_m_max << endl;
			for (int ii=0; ii<=input_m_max; ii++){
				LOG << ii << ": " << bt_cut[ii] << endl;	
			}
		}
//		if (verbose >= 7) LOG << "At acc = " << accuracy[i] << "; i = " << i << endl;
	}
}

int Get_Default_M (double x, int err){
	int e2 = err; int found = 0;
	for (vint_s i=0; i<accuracy.size(); i++){
		if (accuracy[i] >= err)	{
			e2 = accuracy[i];
			found = 1;
			break;
		}
	}
	if (found == 0){
		e2 = accuracy[accuracy.size()-1];
		if (warning>=1) cerr << "Warning! Exceeds the default max relative precision of I0 = 1e-80." << endl;
	}
	
	if (verbose >= 5) LOG << "Getting default M at precision = " << e2   << " for x = " << x << endl;
	if (x == 0){return 0;}
	int start = input_bessel_n;
	if (start < 0) start = 0;
	for (int i = start; i<=input_m_max; i++){
		if (x < bessel_taylor_m_cut[e2][i]){
			return i;	
		}	
	}
	if (warning >=1) cerr << "Warning! The max Taylor-I0 degree is chosen. Result may be inaccurate." << endl;
	return 	input_m_max;
}

void MinDistance(vector<double>& v,  vector<pair <int, int> >& hier){
	if (verbose >= 5) LOG << "Doing hierarchical clustering" << endl;
	int n = v.size();
	vector < pair< pair<int,int>, double> > dist;
	for (int i=0; i<n; i++){
		for (int j=i+1; j<n; j++){
			pair<int,int> id = make_pair(i,j);
			pair<pair<int,int>, double > d0 = make_pair(id, abs(v[i]-v[j]));
			dist.push_back(d0);
		}
	}
	sort(dist.begin(), dist.end(), CMP_PairDB);
	map <int, int> node;
	for (int i=0; i<n; i++){node[i] = i;}
	int con = n;
	for (vp1_s i=0; i<dist.size(); i++){
		pair<int,int> p = dist[i].first;
		int p1 = p.first;
		int p2 = p.second;
		if (node[p1] == node[p2]){continue;}
		pair<int, int> h = make_pair( node[p1], node[p2] );
		hier.push_back(h);
		int a = node[p1];
		int b = node[p2];
		for (map<int,int>::iterator it = node.begin(); it!=node.end();it++){
			if (it->second == a){it->second = con;}	
			if (it->second == b){it->second = con;}	
		}		
		if (con == 2*n - 2){
			break;	
		}
		con ++;
	}
	return;
}

void Exponential (mpf_class expon, mpf_class& e, double error){ // please use positive x
	int sign = 1;
	if (expon < 0){
		sign = -1;
		expon *= (-1);
	}
	
	int e_log2_d = 0;
	mpf_class x(expon);
	if (expon > 1){
		frexp(expon.get_d(), &e_log2_d);
		mpf_div_2exp(x.get_mpf_t(), expon.get_mpf_t(), e_log2_d );
	}

	double est_c;
	int p2 = 1 << e_log2_d;
	if (sign == 1){
		est_c = exp(expon.get_d()) * p2;
	}else{
		est_c = exp( - expon.get_d()) * p2;
	}

	mpf_set_ui(e.get_mpf_t(), 1);
	int k = 1;
	mpf_class inc(1);
	while (inc*est_c > error) {
		inc *= (x/k);
		e += inc;
		k++;
		if (k > fact_max){
			warn_signal=1;
			if (warning == 1) cerr << "Problem at computing exponential. Exceeds FACT_MAX." << endl;
			break;
		}
//		if (verbose >=9) LOG << "k = " << k << "; e = " << e << "; last = " << last << endl;
	}
//	LOG << "e_log2_d = " << e_log2_d << "\texp(x) = " << e << endl;
	for (int i=1; i<= e_log2_d; i++){e = e*e;}	
	if (sign == -1){
		e = 1/e;	
	}

	if (verbose >= 6) LOG << "Exponential of " << expon << " = " << e << "; stop at " << k << endl;
}

void Compute_Next_Error (mpf_class a, mpf_class b, mpf_class c, mpf_class z, int m, mpf_class& e){
	b = b/2.0;
	mpf_class tmp (0);
	mpf_pow_ui(tmp.get_mpf_t(), b.get_mpf_t(), 2*m + 2);
	c = c*tmp; 
	c = c*factorial[2*m+2]/(factorial[m+1]*factorial[m+1]);
	mpf_pow_ui(tmp.get_mpf_t(), a.get_mpf_t(), 2*m + 3);
	c = c/tmp;
	mpf_class one (1);
	mpf_class exps (0);
	z = z*a;
	if (verbose >= 7) LOG << "coef in computing delta = " << c << endl;
	for (int i=0; i<=(2*m+2); i++){
		mpf_pow_ui(tmp.get_mpf_t(), z.get_mpf_t(), i);
		exps += (tmp/factorial[i]);
	}
	if (z < default_my_exp){
		if (verbose >= 5) LOG << "Using Taylor series to compute natural power" << endl;
		Exponential(z, tmp);
		e = c * (one - exps/tmp );
	}else{
		if (verbose >= 5) LOG << "Using C++ exp to compute natural power" << endl;
		e = c * (one - exp(-z.get_d())*exps );
	}
}


void AE (mpf_class& x, mpf_class& e){
	mpf_class inc (1);
	e  = inc;
	for (int i=1; i<80; i++){
		inc = inc * x / i;	
		e += inc;
	}	
	return;
}


