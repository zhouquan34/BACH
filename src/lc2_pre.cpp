/*
   Copyright 09/30/2015 Quan Zhou, Yongtao Guan, Baylor College of Medicine
   This file is part of program bach. You can redistribute it or modify it under the terms of the GNU General Public License. See it for details. You should have received a copy of GPL which is named 'COPYING'.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "lc2_pre.h"

using namespace std;

const int default_min_err_b = 1;
const int default_pre = 20;
const double default_warn_err = 0.1;

int not_compute_p = 0;
int better_bound = 1;
int input_fact_max = -1;
int user_precision = 0;
double small_coef = 0;
double small_r = 0;
int cout_digit = 6;
vector <LinChiSquare*> even_chi; // 2 4 6 8 
string in_file, out_file, log_file;
int next_err = 0;
int limit_err = 0;
int max_p = 50;
int input_omit_deg = -1;
double input_x = 10.0;
double x_for_m = 50.0;
vector <double> chi_coef;
vector <double> cmd_xs;

//ifstream IN;
//ofstream OUT;

// Better_bound: 
// P_bound: even_chisq
// optimal_m: even_chisq

void Compute_Even_ChiSq (int p){
	int esize =  (int) even_chi.size();
	if (p < esize){return;}
	int v_copy = verbose; verbose = 0;	
	for (int i=esize; i<=p; i++){
		int n = 2*(i+1);
		vector <double> c;
		for (int j=1; j<=n; j++) c.push_back(1);
		LinChiSquare* ec = new LinChiSquare(c);  
		ec->Compute_Density();	
		even_chi.push_back(ec);
	}
	verbose = v_copy;
	return;
}

int Scale_Coef(vector<double>& c, double& x, double& scale){
//	sort(c.begin(), c.end());	
	double cmax = c[c.size()-1];
	double sc = 1.0/cmax;
	scale = sc;
	int dim = c.size();
	for (int i=0; i<dim; i++){
		c[i] *= sc;	
	}
	x *= sc;
	int small_mark = -1;
	for (int i=dim-1; i>=0; i--){
		if (c[i] < small_coef){
			small_mark = i;
			break;
		}	
	}

	double cs = 0.0;
	for (int i=0; i<dim; i++) cs += c[i];
	int mark = -1;
	
	double ss = 0.0;
	for (int i=0; i<dim-1; i++){
		ss += c[i];
		if (ss < small_r * cs){ mark = i; }
		else {break;}
	}
	
	int max_mark = mark;
	if (small_mark > max_mark) max_mark = small_mark;
	
	if (max_mark >= 0) c.erase(c.begin(), c.begin()+max_mark+1);

	return max_mark+1;
}

int Estimate_P (double c1, double x){
	if (verbose >= 5) LOG << "estimate pv for x = " << x << ", max coef = " << c1 << endl;
	double x0 = x/c1;
	double p_low = ChiSq_One_PV(x0);
	int mag = -log2(p_low)/log2(10) + 1;
	if (verbose >= 5) LOG << "p =  " << p_low << ", magnitude = " << mag << endl;
	return mag;	
}

void Better_Bound (vector<double> coef, mpf_class x, vector<double>& res){	
//	if (res[1] > 0){return;}
	if (verbose >= 3) {LOG << "Computing tight bounds for x = " << x << endl; LOG << "With coef = "; VecPrint(coef);}
//	int v_copy = verbose;
//	verbose = 0;	
	LinChiSquare* low = NULL;
	LinChiSquare* up = NULL;
	Upper_P(coef, up);
	Lower_P(coef, low);
	res[0] = low->FDF(x);
	res[1] = up->FDF(x);
	delete (low);
	delete (up);
//	verbose = v_copy;
	return;
}

void Upper_P (vector<double>& coef, LinChiSquare*& chi){
	if (verbose >= 5) LOG << "In Upper P" << endl;
	int v_copy = verbose; verbose = 0;
	vector <double> c = coef;
	int cn = c.size();
	if (cn % 2 == 1){ c.insert(c.begin(), c[0]); }
	for (int i=0; i<cn; i+=2){c[i] = c[i+1];}
	chi = new LinChiSquare(c);	
	chi->Compute_Density();	
	verbose = v_copy;
}

void Lower_P (vector<double>& coef, LinChiSquare*& chi){
	if (verbose >= 5) LOG << "In Lower P" << endl;
	int v_copy = verbose; verbose = 0;
	vector <double> c = coef;
	int cn = c.size();
	if (cn % 2 == 1){ c.erase(c.begin(), c.begin()+1); }
	for (int i=1; i<cn; i+=2){c[i] = c[i-1];}
	chi = new LinChiSquare(c);	
	chi->Compute_Density();	
	verbose = v_copy;
}

int P_Bound(vector<double> coef, mpf_class x, vector<double>& res){  // res save estimates; coef already sorted. Even.
	if (res[1] > 0){ 
		if (verbose >= 3) cout << "Pvalue bound computed more than once!" << endl;	
	}
	int cn = coef.size();
	double p_low = ChiSq_One_PV(x.get_d()/coef[cn-1]);
	int tmp_mag = -log10(p_low) + 1;
	if (user_precision>0) default_pv_prec = user_precision + tmp_mag;
	if (verbose >= 3) {LOG << "Estimating x = " << x << " for coef = "; VecPrint(coef);}
	int p = (cn+1)/2 - 1;
	mpf_class x0 (x);
	x0 = x/coef[cn-1];
	res[1] = even_chi[p]->FDF(x0);
	int k = 1;
	for (vdb_s i=2; i<= coef.size(); i+=2){
		double c = coef[ cn - i];	
		x0 = x/c;
		double p_tmp;
		if (x0 > 250+ cn*4 ){  // about pv = 1e-50
			p_tmp = 0;
		}else{
			p_tmp = even_chi[i/2-1] ->FDF(x0);
		}
		if (p_tmp > res[1]) continue;
		if (p_tmp >p_low) { p_low = p_tmp; k = i;}
	}
	res[0] = p_low;
	if (verbose >= 2) LOG << "Lower bound achieved at chi-" << k << endl;		
	
	return -(log10(res[1]*res[0]))/2 + 1;
}

void Optimal_M(vector<double>& coef, mpf_class& x, const int& mpv, vector<int>& best_m){  // res save estimates; coef already sorted. Even.	
	int my_eb =  user_precision + mpv;  // this is log10(gamma) 
	
	vector <double> my_c (coef);	
	if (my_c.size() % 2 == 1) my_c.erase(my_c.begin());
	int cn = my_c.size();
	if (cn < 2) return; 

	// rule delta > gamma/(p*F/2)	
	double x1 = x.get_d();
	best_m.resize(cn/2, 0);
	
	for (int i = 0; i < cn; i += 2){
		mpf_class x2; 
		x2 = x/coef[i];
		double p_tmp = even_chi[0] -> CDF(x2);
		double adjust = cn*p_tmp/2.0;
		int bound = my_eb + log10(adjust);
		if (verbose >= 8) LOG << "i = " << i << "; bound = " << bound << endl;
		if (bound < default_min_err_b) bound = default_min_err_b;
		double b = (my_c[i+1]-my_c[i])/(4.0*my_c[i]*my_c[i+1]);
		best_m[i/2] = Get_Default_M(b*x1, bound);
	}
	return;
}



void Decide_M ( LinChiSquare* mychi, mpf_class& mpf_x, vector<double>& coef, vector<double>& est) {
	double x = mpf_x.get_d();
	int set_m = 0;
	int p_mag = P_Bound(coef, mpf_x, est);
	if (input_bessel_n >= 0){ 
		mychi -> Set_M(input_bessel_n); 
		set_m=1;
		if (verbose >= 2) LOG << "-b is used for choosing M" << endl;
	}
	if (input_ms.size()>0){ 
		mychi -> Set_M(input_ms); 
		set_m=1;
		if (verbose >= 2) LOG << "-B is used for choosing M" << endl;
	}
	if (user_precision >0) { 
		vector<int> my_m;
		Optimal_M(coef, mpf_x, p_mag, my_m);
		mychi->Set_M(my_m); 
		default_pv_prec = user_precision + p_mag;
		set_m=1;
		if (verbose >= 2) LOG << "-e is used for choosing M" << endl;
	}
	if (input_err >0) { 
//		int bound = input_err + log10(coef.size()/2.0) + 0.7; // 0.7 is for correct rounding
		mychi->Use_Default_M(x, input_err); 
		default_pv_prec = input_err;
		set_m=1;
		if (verbose >= 2) LOG << "-E is used for choosing M" <<  endl;
	}
	if (set_m == 0){ 
		mychi->Use_Default_M(x, default_pre);
		if (verbose >= 2) LOG << "Setting M by default"  << endl;
	}

	// settinging the omitting parameter
	if (x<= 1) omit_x_deg = 1;
	else {
		int x_om = log(x)/log(10) + 1; 
		omit_x_deg = x_om;
	}
}

void Start (void){
	if (verbose >= 1){
		log_file = "bach_";
		time_t now; time(&now);
		stringstream tt; tt << now; 
		string s; tt >> s;
		log_file.append(s);
		log_file.append(".log");
		LOG.open(log_file.c_str(), ofstream::out);
	}

	if (input_fact_max != -1) fact_max = input_fact_max;
	if (input_bessel_n >=0 && verbose >=2){LOG << "Input I0 Taylor expansion degree = " << input_bessel_n << endl;}
	if (input_m_max*2 + 2> fact_max){fact_max = input_m_max*2 + 2; }
	if (verbose >=2) LOG << "Precision of floating-point number = " << float_precision << endl;
	if (verbose >= 3){
		LOG << "Max I0 Taylor expansion degree = " << input_m_max << endl;
		LOG << "Verbose mode turned on. Degree = " << verbose << endl;
		if (warning >=1) {LOG << "Warning turned on." << endl; }
		else {LOG << "Warning turned off." << endl;}
	}
	if (not_compute_p > 0){
		if (verbose >= 3) LOG << "P-value and its error bounds will not be computed" << endl;
		pv_bound = 2;
		next_err = 0;
		err_bound = 0;
		limit_err = 0;
	}
	cout.precision(cout_digit);
	LOG.precision(cout_digit);
	mpf_set_default_prec(float_precision);
	Precompute_Fact();
	Compute_Default_M(); // if input_err > 0, changes happen here.
	
//	if (verbose>=2) cout << "Finished precomputation" << endl;	
}

void End(double used_t, int wt ){
	if (wt == 1){
		cerr << "Total used time = " ;
		if (used_t > 2000){cerr << used_t/60.0 << " m" << endl;}
		else{cerr << used_t << " s" << endl;}
	}
	for (int i=0; i< (int) even_chi.size();i++) delete even_chi[i];
}

int Read_Para(int argc, char ** argv){
	int xs_start = 0;
	for (int i=1; i<argc; i++){
//		cout << "Parsing argument " << i << "\t" << argv[i]<< endl;
		if (argv[i][0] == '-'){
			char t_argv = argv[i][1];
			switch(t_argv){							
				case 'a':{
					better_bound = 1; 
					pv_bound = atoi(argv[i+1]);break;
				}
//				case 'b':
//					input_bessel_n = atoi(argv[i+1]);break;
//				case 'B':{
//					string my_m = argv[i+1];
//					replace(my_m.begin(), my_m.end(), ',', ' ');
//					stringstream ms (my_m); int m;
//					while (ms >> m) input_ms.push_back(m);
//					break;
//				}
				case 'c':
					cout_digit = atoi(argv[i+1]);break;
//				case 'd':
//					DOUBLE_PRE = atoi(argv[i+1]);break;
				case 'e':
					user_precision = atoi(argv[i+1]);break;
				case 'E':
					input_err = atoi(argv[i+1]);break;	
				case 'f':
					float_precision = atoi(argv[i+1]);break;
				case 'h':
					return 2; break;	
				case 'i':
					in_file = argv[i+1];break;
				case 'k':
					numeric_interval = atoi(argv[i+1]);break;				
//				case 'l':
//					err_bound = atoi(argv[i+1]);break;				
				case 'm':
					input_m_max = atoi(argv[i+1]);break;
				case 'n':
					input_fact_max = atoi(argv[i+1]);break;
				case 'o':
					out_file = argv[i+1];break;			
//				case 'q':
//					List_Para();break;
				case 'r':
					small_r = atof(argv[i+1]);break;	
				case 's':
					small_coef = atof(argv[i+1]);break;
//				case 'v':
//					verbose = atoi(argv[i+1]);break;
				case 'w':
					warning = atoi(argv[i+1]);break;		
//				case 'z':
//					limit_err = atoi(argv[i+1]);break;
				default:
					return 1;
			}
			i ++;
		}else{
// There should be no whitespace between the comma and the number before it // 
			string x = argv[i];
			if (xs_start == 1 || x.find(",")!=string::npos){
/*				int comma = 0;
				while ( (comma = x.find(",", comma) ) != string::npos ){
					x.replace(comma, 1, " ");
				}
				stringstream xs (x);
				double x1;
				while (xs >> x1){
					cmd_xs.push_back(x1);	
				}
*/
				xs_start = 1;
				stringstream xs(x);
				string x1;
				while (getline(xs, x1, ',')){
					double x0 = atof(x1.c_str());
					cmd_xs.push_back(x0);
//					cout << "X: " << x0 << endl;
				}
			}else{
//				cout << "Coefficients fetched: " << argv[i] << endl;
				chi_coef.push_back(atof(argv[i]));
			}
		}
	}
	
//	cout << "Finished fetching arguments" << endl;
//	cout << chi_coef.size() << endl;
//	cout << cmd_xs.size() << endl;
	if (chi_coef.size() > 0 ){
		if (cmd_xs.size() == 0)	{
			if (chi_coef.size() == 1){
				cerr << "Please input both the cofficients and X" << endl;
				return 1;
			}
			cmd_xs.push_back(chi_coef[chi_coef.size()-1]);	
			chi_coef.erase(chi_coef.end() - 1);
		}
	}
//	cout << chi_coef.size() << endl;
//	cout << cmd_xs.size() << endl;
	return 0;
}


void List_Para(void) {
	cout << "-a: pv_bound" << endl;
	cout << "-b: input_bessel_n" << endl;
	cout << "-B: input_ms" << endl;
	cout << "-c: cout_digit" << endl;
	cout << "-d: DOUBLE_PRE" << endl;
	cout << "-e: user_precision" << endl;
	cout << "-E: input_err" << endl;
	cout << "-f: float_precision" << endl;
	cout << "-h: HELP" << endl;
	cout << "-i: in_file" << endl;
	cout << "-k: numerical_interval" << endl;
	cout << "-l: err_bound" << endl;
	cout << "-m: input_m_max" << endl;
	cout << "-n: input_fact_max" << endl;
	cout << "-o: out_file" << endl;
	cout << "-q: List_Para()" << endl;
	cout << "-r: small_r" << endl;
	cout << "-s: small_coef" << endl;
	cout << "-v: verbose" << endl;
	cout << "-w: warning" << endl;
	cout << "-z: limit_err" << endl;
	exit(EXIT_FAILURE);
}


