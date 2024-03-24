/*
   Copyright 12/02/2015 Quan Zhou, Yongtao Guan, Baylor College of Medicine
   This file is part of program bach. You can redistribute it or modify it under the terms of the GNU General Public License. See it for details. You should have received a copy of GPL which is named 'COPYING'.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "lc2_pre.h"

using namespace std;

void Exit(void){
	cout << "Examples:\n\t./bach 1 0.9 0.8 0.7 10,20\n\t./bach -i test_in.txt\n" << endl;
	cout << "Use -h to read the help file." << endl;
	exit(EXIT_FAILURE);
}

void Help_and_Exit(void){
	cout << "You can input coeffcients and test stat value by command-line arguments or provide an input file to invoke the batch mode." << endl;
	cout << "Please refer to README or User_Manual.pdf for more instructions." << endl;
	cout << endl;
	cout << "Argument [type=default_value]:" << endl;
	cout << "-i [string=\'\']:  (required for batch mode) Input filename." << endl;
	cout << "-o [string=\'\']:  Output filename." << endl;
//	cout << "-v [int=0]:  Verbose mode. If v=1/2/3, write a log file. " << endl;
	cout << "-c [int=6]:  Precision (number of significant digits) for output." << endl;
	cout << "-e [int=5]:  The desired relative error bound (default = 0.001% for batch mode)." << endl;
//	cout << "-a [int=1]:  If a=2, always output the bounds of p-values. If a=0, always not." << endl;
//	cout << "-l [int=1]:  If l=1, output the error bound." << endl;
//	cout << "-s [double=0]:  Omit the coefficients smaller than s when the largest is scaled to 1." << endl;
//	cout << "-r [double=1e-2]:  Omit the small coefs whose sum is smaller than r*(sum of all coefs)." << endl;
//	cout << "-w [int=0]:  If w=1, turn on warning." << endl;
//	cout << "-E [int=0]: The desired absolute error bound." << endl;	
//	cout << "-z [int=0]: z=1 will output the limiting error." << endl;
//	cout << "\nOther precision parameters:" << endl;
//	cout << "-d [int=4]: Effective digits of the input file." << endl;
//	cout << "-k [int=100]: The number of intervals in the numerical integration for an odd number of coefs." << endl;	
//	cout << "-f [int=256]: Set the precision in bits of floating numbers." << endl;
//	cout << "-m [int=80]: The maximum Taylor expansion order of I0." << endl;
//	cout << "-n [int=500]: The maximum factorial allowed, which is also used as the max degree retained in computing density." << endl;
	exit(EXIT_FAILURE);
}

void Output_Coef (vector<double>& c, int o, string& out){
	stringstream coef_str;
	coef_str.precision(DOUBLE_PRE);
	int n = c.size();
	if (n == 1){
		coef_str << c[0];
		coef_str >> out;
		return;
	}
	
	if ( (n-o) % 2 == 1){
		coef_str << c[o] << " "; // may attach a *	
		o = o + 1;	
	}
	for (int i=o; i<( n - 1); i++){
		coef_str << c[i] << " ";	
	}
	coef_str << c[n-1];
//	coef_str >> out;
	out = coef_str.str();
}

void lc2_single (void){
	// Setting default precision //
	if (input_bessel_n <0 && input_ms.size()<1 && input_err <= 0 && user_precision <= 0 ) input_err = 20;  
	Start();
	
	// Compute Even Chisq for Pval bounds//
	Compute_Even_ChiSq( chi_coef.size()/2  );
	if (input_omit_deg > 0) {omit_x_deg = input_omit_deg;}

	sort(chi_coef.begin(), chi_coef.end());	
	double max_x = Max_DB(cmd_xs);
	mpf_class my_max;   Double2MPF(max_x, my_max);
	vector<double> pb(2, 0.0);
	
	LinChiSquare* mychi = new LinChiSquare(chi_coef);	
	Decide_M(mychi, my_max, chi_coef, pb);	
	mychi->Compute_Density();

	for (int i=0; i < (int) cmd_xs.size(); i++){
		mpf_class my_x; 
		Double2MPF (cmd_xs[i], my_x);
		double pv = mychi->FDF(my_x); 
		if (verbose >= 1) {
			LOG << mychi->Output(1);
			double e2 = 0.0;
			LOG << "Pvalue of " << cmd_xs[i] << " = " << pv << endl;
			LOG << "\tLimiting error = " << mychi->Max_Error() << endl;
			e2 = mychi->Compute_Err_Bound(my_x);
			LOG << "\tError bound = " << e2 << endl;	
			vector<double> pb2 (2, 0.0);
			int v_copy = verbose;
			verbose = 0;
			Better_Bound(chi_coef, my_x, pb2); 
			if (pb2[0] > pb[0]){pb[0] = pb2[0];}
			pb[1] = pb2[1];
			verbose = v_copy;
			LOG << "\tLower P bound  = " << pb[0] << "; Upper P bound = " << pb[1] << endl;
		}
		cout << pv << "\t";
	}	
	cout << endl;
	delete(mychi);
	return;
}

void lc2_batch (void){
	if (input_bessel_n <0 && input_ms.size()<1 && input_err <= 0 && user_precision <= 0 ) user_precision = 5;  
	Start();
	
	ifstream IN(in_file.c_str() ); 
	ofstream OUT;
	if (! out_file.empty()){
		OUT.open(out_file.c_str(), ofstream::out);	
		OUT.precision(cout_digit);
//		OUT << "p-val\trow\tomitted\tTaylor-I0\tlastCoef";
//		if (limit_err >=1 ) OUT << "\tlimitingE" ;
//		if (err_bound >=1) OUT << "\terrBound" ;
//		if (pv_bound >=1) OUT << "\tlowerP\tupperP"; 
//		OUT << endl;
		OUT << "p-value\terror-bound\tcoeffcients\tstatistic" << endl;
	}
	string s;
	int row = 0;

	int print_p = 0;
	if (verbose >= 1 || out_file.empty()){print_p = 1;}
	
	while (getline(IN, s)){
		row ++;
		
		/* ------------------   Reading coefficients and X  ---------------------------    */
		vector<double> coef;
		vector<double> xs;
		
		stringstream ss(s); string col; int xs_start = 0;
		while (ss >> col){
			if (xs_start == 1 || col.find(",")!=string::npos){
				xs_start = 1;
				stringstream xx(col); string x1;
				while (getline(xx, x1, ',')){
					double x0 = atof(x1.c_str());
					xs.push_back(x0);
					if (x0 <= 0) {cerr << "Negative value! Exit!" << endl; exit(EXIT_FAILURE);}
				}
			}else{
				double c = atof(col.c_str());
				coef.push_back(c);
				if (c < 0) {cerr << "Negative coefficient! Exit!" << endl; exit(EXIT_FAILURE);}
			}
		}	
		if (xs_start == 0 ){
			if (coef.size() == 1){
				if (warning >= 1 )cerr << "No coefficients specified!" << endl;
				continue;	
			}
			xs.push_back(coef[coef.size()-1]);	
			coef.erase(coef.end() - 1);
		}
		double max_x = Max_DB(xs);
		sort(coef.begin(), coef.end());	 // Important!
		
		/*----------------------- Scale/discard coefficients ---------------------------*/ 
		Compute_Even_ChiSq(coef.size()/2 );
		vector<double> coef_copy = coef;
		double scale;
	
		int omit = Scale_Coef(coef, max_x, scale); // coef, x may change here.
		if (omit != 0) { if (warning>=1) cerr << "Warning! " << omit << " coefficients are omitted" << endl;}
	
		string coef_out; 
		Output_Coef(coef_copy, omit, coef_out);

		mpf_class mpf_max_x;
		Double2MPF(max_x, mpf_max_x, 8);
	
		/*------------------------------ Compute FDF -------------------------*/
		LinChiSquare* mychi = NULL;
		LinChiSquare* low_bound = NULL;
		LinChiSquare* up_bound = NULL;	
		mychi = new LinChiSquare(coef);	
		vector<double> tmp(2, 0.0);
		Decide_M(mychi, mpf_max_x, coef, tmp);
		mychi->Compute_Density();		
		double e1 = mychi->Max_Error();
		
		/*------------------------------ Compute P-values -------------------------*/
		int calculated = 0;
		for (int i=0; i<(int) xs.size(); i++ )	{
			vector<double> pb(2, 0.0);
			vector<double> pb2(2, 0.0);
			
			double x_copy = xs[i]; double x = x_copy;
			mpf_class mpf_x_copy(x_copy);
			mpf_class mpf_x;
			x *= scale;
			if (x<= 1) {omit_x_deg = 1;}
			else {
				int x_om = log(x)/log(10) + 1; 
				omit_x_deg = x_om;
			}
			Double2MPF(x, mpf_x, 8);
				
			if (user_precision > 0) {
				int p_mag = P_Bound(coef, mpf_x, pb);
				default_pv_prec = user_precision + p_mag;
			}

			double pv = mychi->FDF(mpf_x);
			
			/*----------------------- More details ------------------------*/

			double e2 = mychi->Compute_Err_Bound(mpf_x);
			
			int comp = mychi->WarnM();
			if (warn_signal == 1) comp = 1;
			if (e1 > pv || e1 >0.1 ) comp = 1;
			if (e2/pv > default_warn_err) comp = 1;
			
			int pass = 1; 
			double e3 = 1.0;
			if (comp == 1 ){
				if (calculated == 0){
					Upper_P(coef_copy, up_bound);				
					Lower_P(coef_copy, low_bound);		
					calculated = 1;
				}
				pb2[0] = low_bound->FDF(mpf_x_copy);
				pb2[1] = up_bound->FDF(mpf_x_copy);
				if (pb2[0]<pb[0]) pb2[0] = pb[0];
				if ( !(pv>=pb2[0] && pv<= pb2[1] )  ) {	
					pv  = (pb2[0] + pb2[1])/2;
//					pv = pb2[1];
					pass = 0;
					e3 = pv - pb2[0];
				}else{
					e3 = pv - pb2[0];
					if (pb2[1] - pv > e3){
						e3 = pb2[1] - pv;	
					}
				}
			}
			
			/***************  stdout ***************/
			if (print_p >=1){
				if (i != (int) xs.size()-1){
					cout << pv << "\t";
				}else{
					cout << pv << endl;	
				}
			}		

			/************** OUT *************/
			if (out_file.empty()){continue;}
			if (e3 < e2) e2 = e3;
			if (pass == 0){
//				OUT << pv << "\t" << e2 << "\t" << coef_out << "\t" << x_copy << "*" << endl;
				OUT << pv << "\t" << e2 << "\t" << coef_out << "\t" << x_copy  << endl;
			}else{
				OUT << pv << "\t" << e2 << "\t" << coef_out << "\t" << x_copy << endl;
			}
		
		}
		delete(mychi);
		if (calculated == 1){
			delete(low_bound);
			delete(up_bound);
		}
	}
	IN.close();
	if (!out_file.empty()){
		OUT.close();
	}
	return;	
}

int main (int argc, char ** argv){
	if (argc<2){Exit();}  
	warning = 0; verbose = 0; numeric_interval = 0;  // Default for bach
	
	int status = Read_Para(argc, argv);
	if (status == 2) {Help_and_Exit(); }
	if (status == 1) {Exit(); }

	clock_t begin_t = clock();		


	int write_t  = 0;
	if (chi_coef.size() > 0){		
		lc2_single();
		if (! in_file.empty()) cerr << "To invoke the batch mode, please do not input coefficients by command line" <<endl;
	}else{	
		if ( in_file.empty() ) {Help_and_Exit();}
		else { lc2_batch(); }
		write_t = 1;
	}
	
	clock_t end_t = clock();
	double used_t =( (double) end_t - begin_t)/CLOCKS_PER_SEC;
	End(used_t, 0);
	return 0;	
}





