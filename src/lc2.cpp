/*
   Copyright 09/30/2015 Quan Zhou, Yongtao Guan, Baylor College of Medicine
   This file is part of program bach. You can redistribute it or modify it under the terms of the GNU General Public License. See it for details. You should have received a copy of GPL which is named 'COPYING'.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "lc2.h"

using namespace std;

int convolution_k = 0;
int deduct_signal = 0;
int pv_bound = 1;
int float_precision = 256;
int numeric_interval = 0; 
int DOUBLE_PRE = 12;
int omit_on = 1;
int err_bound = 1; 
mpf_class empty_mpf(0);
mpf_class* default_mpf = &empty_mpf;

void Double2MPF (const double& x, mpf_class& f, int prec){
	stringstream ss;
	ss.precision(prec);
	ss << x;
	string ssc(ss.str());
	mpf_set_str(f.get_mpf_t(), ssc.c_str(), 10);
}

LinChiSquare::LinChiSquare(vector<double> coef){	
	odd = 0;
	o = 0.0;
	mpf_set_ui(odd_c.get_mpf_t(), 0);
	coef_sum = 0.0;	
	for (vint_s i=0; i<coef.size(); i++) {coef_sum += coef[i];}
// the sorting of coef should be conducted in the main or elsewhere
	if (coef.size() % 2 == 1){
		odd = 1;
		Double2MPF(coef[0], odd_c);
		coef.erase(coef.begin());
		o = odd_c.get_d();
	}
	if (verbose >=1) {LOG << "Coefficients: "; VecPrint(coef);}	
	dim = 0; simple = 0;
	mpf_set_ui(C.get_mpf_t(), 1);
	int csize = coef.size();
	for (int i=0; i< csize-1; i+=2){
		mpf_class ai(0);
		mpf_class bi(0);
		Double2MPF(coef[i], ai);
		Double2MPF(coef[i+1], bi);
		
		if (verbose >= 4) {
//			LOG << "ai = " ;
//			mpf_out_str( LOG, 10, 20 , ai.get_mpf_t());
//			LOG << "; bi = ";
//			mpf_out_str( LOG, 10, 20 , bi.get_mpf_t());
//			LOG << endl;
			LOG << "ai = " << ai <<  "; bi = " << bi << endl;
		}
		
		dim ++;
		mpf_class ci (4.0*ai*bi);
		mpf_class big_a ( (ai + bi)/ci);
		mpf_class big_b ( (bi - ai)/ci);
		A.push_back( big_a );
		B.push_back( big_b );
		mpf_sqrt(ci.get_mpf_t(), ci.get_mpf_t());
		ci = 1/ci;
		C = C*ci;
		Ci.push_back( ci );
		if (verbose >=2) LOG << "i=" << i << "; Ai = " << big_a << ", Bi = " << big_b << ", Ci = " << ci << endl;
	}	
	c_mag = log10(C.get_d());
	
	if (verbose>=2) LOG << "Scaling constant = " << C << endl;
	pdf = new PolyUniExp();
	cdf = new PolyUniExp();
	fdf = new PolyUniExp();
	
	if (A.size() == 0){
		simple = 1;
		return;
	}
	M.resize(dim,0);
}

LinChiSquare::~LinChiSquare(){
	delete(pdf);	
	delete(cdf);
	delete(fdf);
	for (int i=0; i< (int) BiChi.size(); i++){
		delete (BiChi[i]);
	}
}

double LinChiSquare::Compute_Err_Bound (mpf_class& x){
	if (verbose >= 3) LOG << "ErrBound: Delta for each bi-chi = ";
	mpf_class xe (1); mpf_class fi(0);
	mpf_class error;  mpf_class delta (0); mpf_class tx;
	int tmp_mag = c_mag;
	if (tmp_mag < 0) tmp_mag = 0;
	for (int i=0; i< dim; i++){
		fi = 0;
		BiChi[i]->Compute_Value(x, fi, tmp_mag);
		tx = B[i]*x;
		Compute_I0_Err(tx, M[i], delta);
		error = fi * delta * Ci[i];
		if (verbose >=3) LOG << error << "; ";
		xe *= (error+1);
		if (verbose >= 6) LOG << "\n Finished computing error bound for bi-pdf " << i << endl;
	}	
	if (verbose >= 3 ) LOG << endl;
	xe = xe - 1;
	return (xe.get_d());
}

double LinChiSquare::Compute_Next_Delta (mpf_class& x){
	if (verbose >= 3) LOG << "NextDelta: Error estimates for each bi-chi = ";
	mpf_class xe (1);
	for (int i=0; i< dim; i++){
		mpf_class error;
		Compute_Next_Error(A[i], B[i], Ci[i], x, M[i], error);
		if (verbose >=3) LOG << error << "; ";
		xe *= (error+1);
	}	
	if (verbose >= 3 ) LOG << endl;
	xe = xe - 1;
	return (xe.get_d());
}

string LinChiSquare::Output(int mode){
	if (verbose >= 5){LOG << "Making output" << endl;}
	stringstream m;
	if (mode == 1) m << "Taylor-I0 degrees = ";
	for (int i=0; i< dim; i++){
		if (M[i] < input_m_max){
			m << M[i] << ",";	
		}else{m << M[i] << "*,";}
	}
	if (mode == 1) m << "\t" << "Orphan coef = " << o << endl;
	else m << "\t" << o << "\t";	
	return m.str();
}

void LinChiSquare::Set_M(int b){
	if (verbose >= 5){LOG << "Setting all M = " << b << endl;}
	if ( b > input_m_max){b = input_m_max;}	
	for (int i=0; i<dim; i++){
		if (B[i] == 0.0){
			M[i] = 0;	
		}else{
			M[i] = b;	
		}
	}	
	return;
}

void LinChiSquare::Set_M(vector<int> b){
	if (verbose >= 5){LOG << "Setting M (vector mode)" << endl;}
	if ( (int)b.size() != dim){
		cerr << "The number of input I0 expansion degrees is wrong!\n";
		exit(EXIT_FAILURE);
	}
//	if ( b > fact_max){b = fact_max;}	
	for (int i=0; i<dim; i++){
		if (B[i] == 0.0){
			M[i] = 0;
			continue;
		}
		if (b[i] > input_m_max) M[i] = input_m_max;
		else M[i] = b[i];	
	}	
	return;
}


void LinChiSquare::Use_Default_M(double x, int acc){
	if (verbose >= 5){LOG << "Using default M for x = " << x << " and accuracy = " << acc << endl;}
	for (int i=0; i<dim; i++){
		if (B[i] == 0.0){
			M[i] = 0;
			continue;
		}
		double tmp = B[i].get_d() * x;
		M[i] = Get_Default_M(tmp, acc);
	}	
	return;
}

void LinChiSquare::Compute_Density(void){
	if (simple == 1) {return;}
	if (verbose>=2) {
		LOG << "Taylor-I0 expansion degree: "; 
		for (int i=0; i< dim; i++) LOG << M[i]*2 << "; ";
		LOG << endl;
	}
	vector <double> double_A;
	for (int i=0; i< dim; i++){
		double_A.push_back(A[i].get_d());	
	}
	vector<pair <int, int> > hier;		
	MinDistance(double_A, hier);
	
	if (verbose >= 5){
		LOG << "Integration order by MinDist hierarchical clustering" << endl;
		for (int i=0; i< (int) hier.size(); i++){
			LOG << i << ": " << hier[i].first << ", " << hier[i].second << endl;
		}		
	}
	
	vector<PolyDuoExp*> pdfs;
	for (int i=0; i<dim; i++){
		PolyDuoExp* single = new PolyDuoExp();		
		Single_PDF(i, single);
		pdfs.push_back(single);
		if (err_bound >= 1){
			class PolyUniExp* bi_pdf = new PolyUniExp();	
			bi_pdf->Convert_PDE(single);
			class PolyUniExp* bi_cdf = new PolyUniExp();
			bi_pdf->Integration(bi_cdf);
			BiChi.push_back(bi_cdf);
			delete(bi_pdf);
		}
		if (verbose >= 4){LOG << "Uni PDF of " << i << endl;single->Print(1);}
	}
	
	for (int i=0; i< (int) hier.size(); i++){
		int v1 = hier[i].first;
		int v2 = hier[i].second;
		if (verbose >= 3) LOG << "i = " << i << "; Convolve the " << v1 << "-th and " << v2 << "-th PDFs" << endl;
		Bi_PDF(pdfs[v1], pdfs[v2]);
		pdfs.push_back(pdfs[v1]);
		if (verbose >= 5) {LOG << pdfs.size()-1 << "-th PDF:\n"; pdfs[v1]->Print(1);}
	}
		
	int last = pdfs.size()-1;
	if (verbose >= 3) pdfs[last]->Print_Degree();
	pdf->Convert_PDE(pdfs[last]);
	if (verbose >= 5) {LOG << "pdf (PUE) = \n"; pdf->Print(1);}
	pdf->Integration(cdf);
	if (verbose >= 5) {LOG << "cdf = \n"; cdf->Print(1);}
	cdf->Complementary(fdf, 1/C);
	if (verbose >= 5) {LOG << "fdf = \n"; fdf->Print(1);}
	
	if (verbose >= 2){
		LOG << "PDF = ";
		pdf->Print(1);
		LOG << "CDF = ";
		cdf->Print(1);
		LOG << "FDF = ";
		fdf->Print(1);
	}

	for (int i=0; i<dim; i++){
		delete(pdfs[i]);
	}
	max_e = abs(Max_Error());
	return;	
}

void LinChiSquare::Single_PDF(int v, PolyDuoExp* p1){
	if (verbose >= 5) LOG << "computing single PDF" << endl;
	MultiPoly* zp = new MultiPoly(0);
	zp->Taylor_I0(M[v], B[v]);	
	DuoPoly* start = new DuoPoly(0);
	start->Convert_MP(zp);
	delete(zp);
	p1->Add_Const_DP(start);
	delete(start);
	p1->Multiply_Exp_X(-A[v]);
	return;	
}

void LinChiSquare::Bi_PDF(PolyDuoExp* p1, PolyDuoExp* p2){
	if (verbose >= 7) LOG << "computing bi-PDF" << endl;
	p2->Binomial_Expand();
	if (verbose >= 7) LOG << "Finished binomial expansion" << endl;
	p1->Multiply_PDE(p2);	
	p1->Integration();
	if (verbose >= 7) LOG << "Finished integration" << endl;
	p1->Combine();	
	if (verbose >= 7) LOG << "Finished combining terms" << endl;
	p1->F2X();	
	if (verbose >= 7) LOG << "Finished F2X" << endl;
	if (omit_on >= 1) p1->Omit();
	return;
}

double LinChiSquare::PDF(mpf_class& x){	
	mpf_class p; 
	pdf->Compute_Value(x, p, c_mag);
	p *= C;
	return p.get_d();
}

double LinChiSquare::CDF(mpf_class& x){	
	mpf_class c;
	cdf->Compute_Value(x, c, c_mag);
	c *= C;
	return c.get_d();
}

void LinChiSquare::FDF_Exp_Order(vector<int>& o){
	fdf->Get_Exp_Order(A, o);		
}

double LinChiSquare::FDF(mpf_class& x){
	warn_signal = 0;
	deduct_signal = 0;
	if (verbose >= 5) LOG << "computing p-value of " << x << endl;
	if (x == 0){
		return 1;
	}
	if (x > 130.0 + 4.0*coef_sum){
		if (verbose>=2)	LOG << "x = " << x  << " is big. Be careful with the result" << endl;	
	}
	mpf_class pv (0.0);
	if (odd == 0){
		if (verbose >= 4) LOG << "Computing P-value under the even mode" << endl;
		// decide whether to include epsilon
		vector<double> v_t;
		v_t = fdf->Compute_Value(x, pv, c_mag, 1);
		if (pv > fdf->Constant()  ) {
			pv = pv + fdf->Constant();	
		}else{
			deduct_signal = 1;	
		}
		pv *= C;

	}else{
		if (verbose >= 4) LOG << "Computing P-value under the odd mode" << endl;
		if (simple == 1){pv = ChiSq_One_PV(x.get_d()/o);}
		else{
			pv = FDF_Odd(x);
		}
	}
	if (pv<0 ){ // i think this is no longer possible
		if (warning >=1) cerr << "Pvalue " << pv  << " is smaller than 0! Coercing to be limiting error." << endl;
		pv = max_e;
		warn_signal = 1;
	}	
	if (pv>1 ){ // i think this is no longer possible
		if (warning >=1) cerr << "Pvalue " << pv  << " is bigger than 1!" << endl;
		warn_signal = 1;
	}

//	(*mpv) = pv;
//	if (verbose >=1) LOG << "Limiting error = " << max_e << endl;
	return pv.get_d();
}

double LinChiSquare::Max_Error(void){
	double c0 = fdf->Constant();
	double err = c0 * ( C.get_d() ) ;
	return err;
}

double LinChiSquare::FDF_Odd(mpf_class& z){
	double z0 = z.get_d();
	double pv = 0.0;

	// P(X+Y>z)  = P(Y>z) + \int_0^z f(Y=y)P(X>z-y) dy
	mpf_class tmp(0.0);
	fdf->Compute_Value(z, tmp, c_mag, 1);
	if (tmp > fdf->Constant()  ) {
		tmp = tmp + fdf->Constant();	
	}	
	tmp = tmp * C;
	double p_y = tmp.get_d();
	pv += p_y;

	if (numeric_interval > 0){
		pv += FDF_Part(z0, numeric_interval, 0.0, z0);
	}else{ 
		// split into ten regions. Target relative err < 0.001
		convolution_k = 0;
		vector <double> part;
		double ss = 0.0;
		for (int j=0; j<10; j++){
			double s = FDF_Part(z0,1,z0*j/10, z0*(j+1)/10);
			part.push_back(s);
			ss += s;
		}
		double tmp_pv = pv + ss;
		for (int j=0; j<10; j++){
			double k1 = ( part[j]/tmp_pv ) * 2000.0;
			if (k1 > 1.0){
				int k2 = k1/10.0;
				int k = (k2+1)*10;
				double s = FDF_Part(z0, k, z0*j/10, z0*(j+1)/10);
//				cout << j << ", " << k << ","<< s << endl;
				pv += s;
				convolution_k += k;
			}else{
				pv += part[j];	
				convolution_k += 1;
			}
		}
	}
	return 	pv;
}

double LinChiSquare::FDF_Part(double z, int k, double z1, double z2){
	double n = (double) k;
	mpf_class s1 (0);	
	for (int i=0; i<k; i++){
		double y = (z2-z1) * (i+0.5)/n + z1;
		mpf_class y1 (y);
		mpf_class fy (PDF(y1));
		mpf_class px (ChiSq_One_PV((z-y)/o));
		s1 += fy*px;
	}
	double s = s1.get_d()*(z2-z1)/n;
	return s;
}


int LinChiSquare::WarnM(void){
	for (int i=0; i<dim; i++){
		if (M[i] == input_m_max){
			return 1;	
		}	
	}	
	return 0;
}

