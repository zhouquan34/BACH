/*
   Copyright 09/30/2015 Quan Zhou, Yongtao Guan, Baylor College of Medicine
   This file is part of program bach. You can redistribute it or modify it under the terms of the GNU General Public License. See it for details. You should have received a copy of GPL which is named 'COPYING'.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "algebra.h"

using namespace std;

const double default_max_cpp_e = 500.0;
const int default_mag = -999999;
const int default_omit_mag_bound = 20;
const int default_omit_min = 5;
const double default_exp_zero = 1e-20; // Note that since weights have <16 effective digits (probably only 2-6 are used). difference of exp coefs cannot be 1e-20 or smaller unless by numerical problems  

int default_pv_prec = 10; // this and current status should be put into lc2.cpp. Now just for convenience

int omit_x_deg = 2;
mpf_class neg1 ("-1"); // probably I worry too much
mpf_class pos1 ("1");

mpf_class ExpCoef (const mpf_class& c){
	mpf_class ret;
	if (abs(c)<default_exp_zero){
		ret = 0;
	}else{
		ret = c;
	}
	return ret;
}

MultiPoly::MultiPoly(int n){
	dmax = n;
	coef.resize(dmax+1,0);
}

void MultiPoly::Copy(const class MultiPoly* mp){
	dmax = mp->dmax;	
	coef = mp->coef;
	return;
}

void MultiPoly::Print(int m, int vi){
	char v = vi;
	int start = 0;
	for (int i=0; i<=dmax; i++){
		if (coef[i] == 0.0 ){continue;}
		if (start == 0){
			start = 1;
		}else{
			if (coef[i]>0){LOG << "+" ;}
		}
		if ( coef[i]!=1.0 || i==0){
			LOG << coef[i];	
		}
		if (i!=0){
			LOG << v ;
			if (i!=1){
				LOG << "^" << i; 
			}
		}
	}
	if (m==1){LOG << endl;}
	return;
}

void MultiPoly::Add_Term(int i, const mpf_class& c){
	if (verbose >= 9) LOG << "MP add term. Degree = " << i << "; Coef = " << c << endl;
	if (i>dmax){
		coef.resize(i+1, 0.0);	
		dmax = i;
	}
	coef[i] += c;
	return;
}

void MultiPoly::Add_MP(const class MultiPoly* mp){
	for (int i=0; i<=mp->dmax; i++){
		Add_Term(i, mp->coef[i]);
	}	
	return;
}

void MultiPoly::Delete_Term(int i){
	if (verbose >= 8) LOG << "MP delete term " << i << endl;
	if (i>dmax){return;}	
	if (i==dmax){
		coef.erase(coef.begin()+dmax);	
		dmax --;
	}else{
		coef[i] = 0;	
	}
	if (verbose >= 8) LOG << "Max degree = " << dmax << endl;
	return;
}

void MultiPoly::Multiply_Scalar(const mpf_class& c){
	if (verbose >= 9) LOG << "MP multiply  " << c << endl;
	for (int i=0; i<=dmax; i++){
		coef[i]	 *= c;
	}
	return;
}

void MultiPoly::Multiply_MP(const class MultiPoly* mp){
	int dsum = dmax + mp->dmax;
	vector<mpf_class> pc(dsum+1, 0.0);	
	for (int i=0; i<=dmax; i++){
		for (int j=0; j<=mp->dmax; j++){
			int d = i + j;
			mpf_class c (coef[i]*mp->coef[j]);
			pc[d] += c;
		}	
	}
	dmax = dsum;
	coef = pc;
	return;
}

void MultiPoly::Convert_DP(const class DuoPoly* dp){
	if (verbose >= 8) LOG << "Converting MP to DP" << endl;
	dmax = 0;
	coef.clear();
	for (int i=0; i<=dp->dmax; i++){
		for (int j=0; j<=dp->coef[i]->dmax; j++){
			Add_Term(i+j,  dp->coef[i]->coef[j]);
		}	
	}
	return;
}

void MultiPoly::Binomial_Expand (class DuoPoly* dp){// expand as (p-x)^n
	if (verbose >= 8) LOG << "Binomial expansion of MP to DP" << endl;
	for (int i=0; i<=dmax; i++){
		if (verbose >= 8) LOG << "i = " << i << endl; 
		mpf_class c (coef[i]);
		if (c == 0.0 ){continue;}
		for (int j=0; j<=i; j++){
			if (verbose >= 8) LOG << "j = " << j << endl; 
			mpz_class coef;
			BinCoef(i, j, coef);
			int sign = 1;
			if (j % 2 == 1) sign = -1;
			coef = coef * sign;
			mpz_class binc(coef);  // right. from mpz to mpf
			int k = i-j;
			MultiPoly* tmp = new MultiPoly(k);
			tmp->Add_Term(k, c*binc); // right. c is mpf, binc is mpz. the product is mpf.
			dp->Add_MP(j, tmp);
			delete(tmp);
		}
	}	
	return;	
}

void MultiPoly::Taylor_I0 (int n, const mpf_class& b){
	if (verbose >= 6) LOG << "Setting MP to be Taylor expansion of I0 with n = " << n << " and coef = " << b << endl;
	if ( b == 0.0 ){
		Add_Term(0, 1);
		return;
	}
	if (n > fact_max){n = fact_max;}
	for (int m=0; m<=n; m++){
		mpf_class nu(b/2.0);
		mpf_pow_ui(nu.get_mpf_t(), nu.get_mpf_t(), 2*m);
		mpf_class k (nu/(factorial[m]*factorial[m]));
		Add_Term(2*m, k);		
	}
	return;	
}

int MultiPoly::If_Zero(void){
	for (int i=0; i<=dmax; i++){
		if ( coef[i] != 0.0){
			return 0;	
		}
	}	
	return 1;
}

int MultiPoly::Magnitude(void){
	int m = default_mag;
	for (int i=0; i<=dmax; i++){
		if ( coef[i] == 0.0){continue;}
		int mi = log(coef[i].get_d())/log(10);
		if (mi>m){m = mi;}
	}
	return m;
}

void MultiPoly::Compute_Value(const mpf_class& x, mpf_class& v){
	for (int i=0; i<=dmax; i++){
		mpf_class xp;
		mpf_pow_ui (xp.get_mpf_t(), x.get_mpf_t(), i);
		v += (xp*coef[i]);	
		if (verbose >= 7) LOG <<  "At the end of degree " << i << ", v = " << v << endl;
	}
}

DuoPoly::DuoPoly(int n){
	dmax = n;	
	for (int i=0; i<=dmax; i++){
		MultiPoly* mpc = new MultiPoly(0);	
		coef.push_back(mpc);
	}
}

DuoPoly::~DuoPoly(){
	for (int i=0; i<=dmax; i++){delete(coef[i]);}	
}

void DuoPoly::Copy(const class DuoPoly* dp){
	for (int i=0; i<=dmax; i++){delete(coef[i]);}
	coef.clear();
	dmax = dp->dmax;	
	for (int i=0; i<=dmax; i++){
		MultiPoly* mpc = new MultiPoly(0);	
		coef.push_back(mpc);
	}
	Add_DP(dp);
	return;
}

void DuoPoly::Print(int m, int vi){
	char v = vi;
	LOG << "[";
	coef[0]->Print(0, 102);
	LOG << "]";
	for (int i=1; i<=dmax; i++){
		LOG << "+[";		
		coef[i]->Print(0, 102);
		LOG << "]" << v << "^" << i;
	}
	if (m==1){LOG << endl;}
	return;
}

void DuoPoly::Print_Dim(void){	
	for (int i=0; i<=dmax; i++){
		LOG << i << ":" << coef[i]->dmax << "; ";	
	}	
	LOG << endl;
	return;
}

void DuoPoly::Add_Term(int i, int j, const mpf_class& c){
	if (i>dmax){
		for (int k=dmax+1; k<=i; k++){
			MultiPoly* mpc = new MultiPoly(0);
			coef.push_back(mpc);
		}
		dmax = i;
	}
	coef[i]->Add_Term(j, c);
	return;
}

void DuoPoly::Add_MP(int i, const class MultiPoly* mp){
	for (int j=0; j<=mp->dmax; j++){
		Add_Term(i, j, mp->coef[j]);
	}
	return;
}

void DuoPoly::Add_DP(const class DuoPoly* dp){
	for (int i=0; i<=dp->dmax; i++){
		Add_MP(i, dp->coef[i]);
	}
	return;
}

void DuoPoly::Multiply_Scalar(const mpf_class& c){
	for (int i=0; i<=dmax; i++){
		coef[i]->Multiply_Scalar(c);
	}
	return;
}

void DuoPoly::Multiply_MP_X(const class MultiPoly* mp){
	int dsum = dmax + mp->dmax;
	vector<class MultiPoly*> pc;
	for (int i=0; i<=dsum; i++){ 
		MultiPoly* mpc = new MultiPoly(0);
		pc.push_back(mpc);
	}
	for (int i=0; i<=dmax; i++){
		for (int j=0; j<=mp->dmax; j++){
			mpf_class c (mp->coef[j]);
			if (c == 0.0){continue;}
			MultiPoly tmp(0);		
			tmp.Copy(coef[i]);
			tmp.Multiply_Scalar(c);
			pc[i+j]->Add_MP(&tmp);
		}	
	}		
	for (int i=0; i<=dmax; i++){delete(coef[i]);}
	coef = pc;
	dmax = dsum;
	return;
}

void DuoPoly::Multiply_DP(const class DuoPoly* dp){
	int dsum = dmax + dp->dmax;
	vector<class MultiPoly*>  pc;
	for (int i=0; i<=dsum; i++){ 
		MultiPoly* mpc = new MultiPoly(0);
		pc.push_back(mpc);
	}
	for (int i=0; i<=dmax; i++){
		for (int j=0; j<=dp->dmax; j++){
			MultiPoly tmp(coef[i]->dmax);		
			tmp.Copy(coef[i]);
			tmp.Multiply_MP(dp->coef[j]);
			pc[i+j]->Add_MP(&tmp);
		}	
	}		
	for (int i=0; i<=dmax; i++){delete(coef[i]);}
	coef = pc;
	dmax = dsum;
	return;	
}

void DuoPoly::Convert_MP(const class MultiPoly* mp){
	for (int i=0; i<=dmax; i++){delete(coef[i]);}
	coef.clear();
	dmax = mp->dmax;	
	for (int i=0; i<=dmax; i++){
		MultiPoly* mpc = new MultiPoly(0);	
		coef.push_back(mpc);
	}	
	for (int i=0; i<=dmax; i++){
		coef[i]->Add_Term(0, mp->coef[i]);
	}
	return;
}

void DuoPoly::Binomial_Expand(void){// expand both X and P as (z-x)
	DuoPoly bi_prod (0); 
	for (int i=0; i<=dmax; i++){
		DuoPoly bpower(0);	
		MultiPoly tmp(i); 
		// perhaps dangerous to write Add_Term(i, 1)
		tmp.Add_Term(i, pos1);
		tmp.Binomial_Expand(&bpower);
		DuoPoly bpoly (0);
		tmp.Copy(coef[i]);
		tmp.Binomial_Expand(&bpoly);
		bpoly.Multiply_DP(&bpower);
		bi_prod.Add_DP(&bpoly);
	}
	Copy(&bi_prod);
	return;
}

void DuoPoly::Gamma_Int (const mpf_class& a, class DuoPoly* cdp, class DuoPoly* edp){ 
	if (verbose >= 6) LOG << "Gamma int of DP with a = " << a  << endl; 
	for (int i=0; i<=dmax; i++){
		MultiPoly tmp (0);
		tmp.Copy(coef[i]);
		if (i>fact_max){break;}
		mpf_class c0 (factorial[i]);
		mpf_class ap ;
		mpf_pow_ui (ap.get_mpf_t(), a.get_mpf_t(), i+1);
		c0 = c0/ap;
		tmp.Multiply_Scalar(c0);
		cdp->Add_MP(0, &tmp);	
	
		for (int j=0; j<=i; j++){
			tmp.Copy(coef[i]);
			mpf_class c1 (-factorial[i]);
			c1 = c1/factorial[j];
			mpf_pow_ui (ap.get_mpf_t(), a.get_mpf_t(), i+1-j);
			if (verbose >= 8) LOG << "i = " << i << "; j = " << j << "; ap = " << ap << endl;
			c1 = c1 /ap;	
			tmp.Multiply_Scalar(c1);
			edp->Add_MP(j, &tmp);	
		}
	}
	return;
}

void DuoPoly::Power_Int (class DuoPoly* dp){	// this is dangerous. dp must be empty.
	if (verbose >= 6) LOG << "Power int of DP" << endl; 
	for (int i=0; i<=dmax; i++){
		dp->Add_MP(i+1, coef[i]);
		mpz_class zi (i+1);
		mpf_class c ("1");
		c = c/zi;
		dp->coef[i+1]->Multiply_Scalar(c);
	}	
	return;	
}

void DuoPoly::Transform_MP(class MultiPoly* mp){
	for (int i=0; i<=dmax; i++){
		for (int j=0; j<=coef[i]->dmax; j++){
			mp->Add_Term(i+j, coef[i]->coef[j]);
		}	
	}
	return;	
}

int DuoPoly::If_Zero(void){
	for (int i=0; i<=dmax; i++){
		if (coef[i]->If_Zero() == 0){
			return 0;
		}	
	}
	return 1;
}

void DuoPoly::F2X(void){
	DuoPoly tmp(0);
	for (int i=0; i<=dmax; i++){
		for (int j=0; j<=coef[i]->dmax; j++){
			tmp.Add_Term(i+j, 0, coef[i]->coef[j]);
		}
	}
	Copy(&tmp);
	return;	
}

void DuoPoly::Omit(void){  // PDF_MIN is not really needed since we only care about big x ... and we set a lower bound on omitting	
	if (verbose >= 5) {
		if (dmax >= fact_max)	LOG << "Omitting " << dmax - fact_max + 1 << " terms that exceed factorial max" << endl;
	}
	for (int i=dmax; i>= fact_max; i--){
		delete(coef[i]);
		dmax -- ;
	}

	int umax = default_mag;
	if (verbose >=6) LOG << "Beginning of Omit, dmax = " << dmax << endl;
	for (int i=0; i<=dmax; i++){
		if (coef[i]->If_Zero() == 1){continue;}
		signed long int e2;
		mpf_get_d_2exp(&e2, coef[i]->coef[0].get_mpf_t());
		int um = ( (double) e2)/log2(10.0) + i * omit_x_deg;
		if (um > umax){ umax = um;}
	}
	if (verbose >= 7) LOG << "magnitude = " << umax  << endl;
	for (int i=dmax; i>= default_omit_min; i--){			
		if (coef[i]->If_Zero() == 1){continue;}				
		signed long int e2;
		mpf_get_d_2exp(&e2, coef[i]->coef[0].get_mpf_t());
		int um = ( (double) e2)/log2(10.0) + i * omit_x_deg;	
		if (verbose >= 8) LOG << "order " << i  << "; um = " << um << "; coef = "<<  coef[i]->coef[0] << endl;
		if ( umax - um > default_omit_mag_bound  ){
			if (verbose >= 7) LOG << "Omit order " << i  << "; um = " << um << "; coef = "<<  coef[i]->coef[0] << endl;
			coef[i]->coef[0] = 0;	
		}else{
			break;	
		}
	}
	for (int i=dmax; i>= default_omit_min; i--){
		if (coef[i]->coef[0] == 0){
			delete(coef[i]);
			dmax -- ;
		}else{
			break;	
		}
	}
	if (verbose >=6) LOG << "End of Omit, dmax = " << dmax << endl;
	return;
}

int DuoPoly::Max_Deg(void){
	return dmax;	
}

DuoExpPoly::DuoExpPoly(){
	mpf_set_ui(pcoef.get_mpf_t(), 0);
	mpf_set_ui(xcoef.get_mpf_t(), 0);
	epoly = new DuoPoly (0);	
}


DuoExpPoly::DuoExpPoly(const mpf_class& cx){
	xcoef = cx;
	mpf_set_ui(pcoef.get_mpf_t(), 0);
	epoly = new DuoPoly (0);	
}

DuoExpPoly::DuoExpPoly(const mpf_class& cx, const mpf_class& px){
	xcoef = cx;
	pcoef = px;
	epoly = new DuoPoly(0);
}

DuoExpPoly::~DuoExpPoly(){
	delete(epoly);
}

void DuoExpPoly::Initialize(void){
	epoly->Add_Term(0, 0, pos1);
	return;
}

void DuoExpPoly::Copy(const class DuoExpPoly* dep){
	xcoef = dep->xcoef;
	pcoef = dep->pcoef;
	epoly->Copy(dep->epoly);
	return;	
}

void DuoExpPoly::Print(int m, int vi){
	char v = vi;
	LOG << "e^(";
	if (pcoef!=0){
		LOG << pcoef << "f+" ; 
	}	
	LOG << xcoef << v << "){";
	epoly->Print();
	LOG << "}";
	if (m==1){LOG << endl;}
	return;	
}

void DuoExpPoly::Add_DP(const class DuoPoly* dp){
	epoly->Add_DP(dp);
	return;	
}

void DuoExpPoly::Multiply_Exp_X(const mpf_class& a){
	xcoef += a;
	return;
}

void DuoExpPoly::Multiply_MP_X(const class MultiPoly* mp){
	epoly->Multiply_MP_X(mp);
	return;	
}

void DuoExpPoly::Binomial_Expand(void){
	epoly->Binomial_Expand();
	if (pcoef!= 0.0){
		if (warning >=1) cerr << "Waring in integration! Exp-poly coef is not zero! " << pcoef << endl;	
	}
	pcoef += xcoef;   // When do binomial expansion, pcoef should always be 0. 	
	xcoef *= neg1;
	return;	
}

void DuoExpPoly::Integration(class PolyDuoExp* dep){
	if (verbose >= 6)  LOG << "DEP integration. pcoef = " << pcoef << "; xcoef = " << xcoef << endl;
	if (pcoef != 0.0){
		if (xcoef == 0.0){
			DuoExpPoly* int_ep = new DuoExpPoly(pcoef);
			epoly->Power_Int(int_ep->epoly);
			dep->Add_DEP(int_ep);
			delete(int_ep);
		}else{
			DuoExpPoly* int_e = new DuoExpPoly(xcoef+pcoef);
			DuoExpPoly* int_c = new DuoExpPoly(pcoef);
			epoly->Gamma_Int(xcoef*neg1, int_c->epoly, int_e->epoly);
			dep->Add_DEP(int_e);
			dep->Add_DEP(int_c);	
			delete(int_e);
			delete(int_c);
		}
	}else{
		if ( xcoef == 0.0){
			if (verbose >= 7) LOG << "DEP Power int" << endl;
			DuoPoly* int_dp = new DuoPoly(0);
			epoly->Power_Int(int_dp);
			dep->Add_Const_DP(int_dp);
			delete(int_dp);
		}else{
			if (verbose >= 7) LOG << "Gamma int" << endl;
			DuoExpPoly* int_e = new DuoExpPoly(xcoef+pcoef);
			DuoPoly* int_c = new DuoPoly(0);
			epoly->Gamma_Int(xcoef*neg1, int_c, int_e->epoly);
			dep->Add_DEP(int_e);
			dep->Add_Const_DP(int_c);	
			delete(int_e);
			delete(int_c);
		}
	}
	return;
}

PolyDuoExp::PolyDuoExp(){
	DuoExpPoly* cep = new DuoExpPoly(0.0);
	dexp.push_back(cep);
	term = 1;
}

PolyDuoExp::~PolyDuoExp(){
	for (int i=0; i<term; i++){
		delete(dexp[i]);
	}
}

void PolyDuoExp::Copy(const class PolyDuoExp* pde){
	if (pde->term<1){
		if (warning >= 1) cerr << "Copying empty PolyDuoExp" << endl; 
		return;
	}	
	for (int i=0; i<term; i++){
		delete(dexp[i]);
	}		
	dexp.clear();
	for (int i=0; i<pde->term; i++){
		DuoExpPoly* ep = new DuoExpPoly();
		ep->Copy(pde->dexp[i]);
		dexp.push_back(ep);					
	}
	term = pde->term;
	return;
}

void PolyDuoExp::Print(int m, int vi){
	dexp[0]->Print(0, vi);	
	for (int i=1; i<term; i++){
		LOG << "\n    +";
		dexp[i]->Print(0, vi);
	}		
	if (m==1){LOG << endl;}
	return;
}

void PolyDuoExp::Add_Term(const mpf_class& ec){
	DuoExpPoly* ep = new DuoExpPoly(ec);
	ep->Initialize();
	dexp.push_back(ep);
	term ++;
	return;
}

void PolyDuoExp::Add_Term(const mpf_class& ec, const mpf_class& pc){
	DuoExpPoly* ep = new DuoExpPoly(ec, pc);
	ep->Initialize();
	dexp.push_back(ep);
	term ++;
	return;
}

void PolyDuoExp::Add_DEP(class DuoExpPoly* dep){
	DuoExpPoly* ep = new DuoExpPoly(0.0);
	ep->Copy(dep);
	dexp.push_back(ep);
	term ++;
	return;	
}

void PolyDuoExp::Add_Const_DP(const class DuoPoly* dp){
	dexp[0]->Add_DP(dp);
	return;
}

void PolyDuoExp::Multiply_Exp_X(const mpf_class& c){
	if (verbose >= 7) LOG << "PDE Multiply_EXP_X start;  term = " << term << endl;
	for (int i=1; i<term; i++){
		dexp[i]->Multiply_Exp_X(c);
	}
	DuoExpPoly* nep = new DuoExpPoly();
	nep->Copy(dexp[0]);
	nep->Multiply_Exp_X(c);
	dexp.push_back(nep);
	term ++;
	
	DuoExpPoly* empty = new DuoExpPoly();
	delete(dexp[0]);	
	dexp[0] = empty;
	if (verbose >= 7) LOG << "PDE Multiply_EXP_X end;  term = " << term << endl;
	return;
}

void PolyDuoExp::Multiply_PDE(const class PolyDuoExp* pde){
	if (verbose >= 7) LOG << "PDE Multiply_PDE start;  term = " << term << endl;
	PolyDuoExp* prod = new PolyDuoExp();	
	for (int i=0; i<term; i++){
		if (dexp[i]->epoly->If_Zero()==1){continue;}
		for (int j=0; j<pde->term; j++){
			if (pde->dexp[j]->epoly->If_Zero()==1){continue;}
			DuoExpPoly* tmp = new DuoExpPoly(ExpCoef(dexp[i]->xcoef + pde->dexp[j]->xcoef), ExpCoef(dexp[i]->pcoef + pde->dexp[j]->pcoef));
			tmp->Initialize();
			tmp->epoly->Multiply_DP(dexp[i]->epoly);
			tmp->epoly->Multiply_DP(pde->dexp[j]->epoly);
			prod->Add_DEP(tmp);
			delete(tmp);
		}	
	}
	Copy(prod);
	delete(prod);
	if (verbose >= 7) LOG << "PDE Multiply_PDE end;  term = " << term << endl;
	return;	
}

void PolyDuoExp::Multiply_MP_X(const class MultiPoly* mp){
	for (int i=0; i<term; i++){
		dexp[i]->Multiply_MP_X(mp);
	}
	return;
}

void PolyDuoExp::Binomial_Expand(void){
	if (verbose >= 7) LOG << "PDE binomial expansion;  term = " << term << endl;
	for (int i=0; i<term; i++){
		if (verbose >= 7) LOG << "Expanding " << i << "-th PDE" << endl;
		if (verbose >= 7) dexp[i]->Print(1);
		dexp[i]->Binomial_Expand();
	}
	if (verbose >= 7) LOG << "PDE binomial expansion end;  term = " << term << endl;
	return;
}	
		
int PolyDuoExp::Combine(void){ // Dangerous! this should be performed only when pcoef=0
	if (verbose >= 7) LOG << "PDE combine start;  term = " << term << endl;
	int cb = 0;	
	vector<int> xcs;
	vector<int> dup;
	for (int i=0; i<term; i++){	
		if (dexp[i]->pcoef != 0.0){
			if (warning >=1) cerr << "Warning in combining! Exp-poly = " << dexp[i]->pcoef << endl;	
		}
		double xc = dexp[i]->xcoef.get_d();
		int found = 0;
		for (vint_s j=0; j<xcs.size(); j++){
//			int mag = log10(abs(xc)) + 1;
//			double diff  = abs(dexp[xcs[j]]->xcoef.get_d() - xc) / pow(10.0,mag);
			double diff  = abs(dexp[xcs[j]]->xcoef.get_d() - xc) ;
			if (verbose >= 9) LOG << "i = " << i << "; j = " << j << "; diff = " << diff << endl;
			if (diff < default_exp_zero){
					if (verbose >= 7) LOG << "Combining i = " << i << "; j = " << j << "; diff = " << diff << endl;
					dexp[xcs[j]]->epoly->Add_DP(dexp[i]->epoly);
					delete(dexp[i]);
					dup.push_back(i);
					cb ++;
					found = 1;
					break;
			}
		}
		if (found == 0){
			xcs.push_back(i);
		}	
	}
	if (verbose >= 6) LOG << "Deleting " << dup.size() << " terms" << endl;
	term -= dup.size();
	for (int i=dup.size()-1; i>=0; i--){
		dexp.erase(dexp.begin()+dup[i]);
	}
	if (verbose >= 7) LOG << "PDE combine end;  term = " << term << endl;
	return cb;
}

void PolyDuoExp::Integration(void){
	if (verbose >= 5) LOG << "PDE integration start;  term = " << term << endl;
	PolyDuoExp* integral = new PolyDuoExp();	
	for (int i=0; i<term; i++){
		if (verbose >= 6) LOG << "Integrating term " << i << endl;
		dexp[i]->Integration(integral);
		if (verbose >= 6) integral->Print(1);
	}
	Copy(integral);
	delete(integral);
	if (verbose >= 7) LOG << "PDE integration end;  term = " << term << endl;
	return;
}

void PolyDuoExp::F2X(void){
	if (verbose >= 7) LOG << "PDE F2X start;  term = " << term << endl;
	for (int i=0; i<term; i++){//pcoef should be 0	
		dexp[i]->epoly->F2X();
	}
	if (verbose >= 7) LOG << "PDE F2X end;  term = " << term << endl;
	return;	
}

void PolyDuoExp::Omit(void){
	if (verbose >= 7) LOG << "PDE Omit start;  term = " << term << endl;
	for (int i=0; i<term; i++){
		if (verbose >=6) LOG << "To Omit DEP " << i << endl;
		dexp[i]->epoly->Omit();	
	}	
	if (verbose >= 7) LOG << "PDE Omit end;  term = " << term << endl;
	return;
}

int PolyDuoExp::Max_Deg(void){
	int maxd = 0;
	for (int i=0; i<term; i++){
		int d = dexp[i]->epoly->Max_Deg();	
		if (d>maxd){maxd = d;}
	}
	return maxd;
}

void PolyDuoExp::Print_Degree(void){
	LOG << "Degree = " ;
	for (int i=0; i<term; i++){
		LOG << dexp[i]->epoly->Max_Deg() << ", ";	
	}
	LOG << endl;
	return;
}

PolyUniExp::PolyUniExp(){
	term = 1;
	ecoef.resize(1,0.0);
	MultiPoly* cpoly = new MultiPoly(0);
	coef.push_back(cpoly);
}

PolyUniExp::~PolyUniExp(){
	for (int i=0; i<term; i++){
		delete(coef[i]);	
	}	
}

void PolyUniExp::Print(int m, int vi){
	char v = vi;
	if (coef[0]->If_Zero() == 0){
		LOG << "[";
		coef[0]->Print();
		LOG << "]\n";
		if (term>1){LOG << "    +";}
	}
	for (int i=1; i<term; i++){
		if (i>1){LOG << "    +";}
		LOG << "e^(" << ecoef[i] << v << ")";
		LOG << "[";
		coef[i]->Print();
		LOG << "]\n";
	}
	return;	
}

void PolyUniExp::Convert_PDE(class PolyDuoExp* pde){
	for (int i=0; i<term; i++){
		delete(coef[i]);	
	}
	term = pde->term;
	ecoef.resize(term, 0.0);
	coef.clear();
	for (int i=0; i<term; i++){
		ecoef[i]=pde->dexp[i]->xcoef; //pcoef should be 0	
		MultiPoly* tmp = new MultiPoly(0);
		pde->dexp[i]->epoly->Transform_MP(tmp);
		coef.push_back(tmp);
	}
	return;	
}

void PolyUniExp::Integration(class PolyUniExp* pue){ // pue must be empty and current pue has no constant term
	if (verbose >= 4) LOG << "PUE integration" << endl;
	for (int i=1; i<term; i++){
		if (verbose >= 8) LOG << "In PUE-integration. i=" << i << endl;
		pue->ecoef.push_back(ecoef[i]);
		pue->term++;
		MultiPoly* cf = new MultiPoly(0);
		pue->coef.push_back(cf);
		mpf_class a (ecoef[i] * neg1);
		for (int j=0; j<=coef[i]->dmax; j++){
			mpf_class c0 (coef[i]->coef[j]);
			if ( c0 == 0.0) {continue;}
			if (j>fact_max){break;}
			mpf_class const_term(c0*factorial[j]);
			mpf_class den;
			mpf_pow_ui(den.get_mpf_t(), a.get_mpf_t(), j+1);
			const_term = const_term/den;
			if (verbose >= 8) {
				long int e2 = 0;
				mpf_get_d_2exp(&e2, const_term.get_mpf_t());
				LOG << "j = " <<  j  << "; const = " << const_term << "; precision = "<<  const_term.get_prec() << "; 2exp = "<<  e2 <<  endl;
			}
			pue->coef[0]->Add_Term(0, const_term);
			int t = (pue->term - 1);
			for (int k=0; k<=j; k++){
				mpf_class c1 (-c0*factorial[j]);
				c1 = c1/factorial[k];  
				mpf_pow_ui(den.get_mpf_t(), a.get_mpf_t(), j+1-k);
				c1 = c1/den;
				pue->coef[t]->Add_Term(k, c1);
			}
		}
	}
	return;
}

void PolyUniExp::Complementary(class PolyUniExp* pue, const mpf_class& c){// pue must be empty
	mpf_class c0 (c - coef[0]->coef[0]);
	if (verbose >= 4) LOG << "In complementary, error = " << c0 << endl;
	pue->coef[0]->Add_Term(0,c0);
	for (int i=1; i<term; i++){
		pue->ecoef.push_back(ecoef[i]);
		pue->term++;
		MultiPoly* cf = new MultiPoly(0);
		cf->Copy(coef[i]);
		cf->Multiply_Scalar(neg1);
		pue->coef.push_back(cf);
	}
	return;
}


vector<double> PolyUniExp::Compute_Value(const mpf_class& x, mpf_class& v, int lc2_c_mag, int start){
	vector<double> ds;
	for (int i=0; i<start; i++) ds.push_back(0.0);
	if (verbose >= 5) LOG << "Computing value of x = " << x  << endl;
	mpf_set_ui(v.get_mpf_t(), 0);
	if (verbose >= 5) LOG << "pv_prec = "  << default_pv_prec << endl;
	
	int call_exp = 0;
	for (int i=start; i<term; i++){
		if (verbose >= 7) LOG << "Compute term " << i <<  endl;
		mpf_class p (0.0);
		coef[i]->Compute_Value(x, p);
		mpf_class e (ecoef[i]*x);
		int cpp_prec =  - (e.get_d()/log(10)) + 15;
		int prec = log10(abs(p.get_d())) + 1 + default_pv_prec + lc2_c_mag ;
		mpf_class delta(0);
		if (e > default_max_cpp_e){
			warn_signal = 1;
			if (warning >=1) cerr << "Warning! Too big value for computing exp!\n" << endl;
			cpp_prec = - 9999;	
		}
//		if (abs(x)>1000) LOG << "i = " << i << "; p = " << p << "; e = " << e << "; default_pv_prec = " << default_pv_prec << "; cpp_prec = " << cpp_prec << "; prec = " << prec << endl;
		if (cpp_prec > prec){
			double cpp_e = exp(e.get_d());
			delta = p*cpp_e;
		}else{
			if (verbose >= 5) LOG << "Using own exp function; cpp_prec = " << cpp_prec << "; prec = " << prec << endl;
			call_exp ++;
			mpf_class ee (0);
			Exponential(e, ee, pow(10.0, -prec));
			delta = p*ee;
//			LOG << ee << endl;
		}
//		mpf_class tmp_ee (0);
//		AE(e, tmp_ee);
//		LOG << tmp_ee << endl;

		v += delta;
		if (verbose >= 4) LOG << "Value of ExpPoly term " << i << "  = " << delta << endl;
		ds.push_back(delta.get_d());		
	}	
	if (verbose >=4 && call_exp>0) LOG << "Calling my own exp functions " << call_exp << " times" << endl;
	return ds;
}


double PolyUniExp::Constant(void){
	return coef[0]->coef[0].get_d();	
}

void PolyUniExp::Get_Exp_Order(vector<mpf_class>& a, vector<int>& o){
	int n = a.size();
	o.resize(n, 0);
	for (int i=0; i<n; i++){
		for (int j=1; j<=n; j++){
//			if (a[i]*neg1 == ecoef[j]){
			if (abs(a[i]*neg1 - ecoef[j])/abs(ecoef[j])< default_exp_zero) {
				o[i] = j;
				break;
			}
		}
		if (o[i] == 0){
			if (warning >= 1) cerr << "Wrong! a_" << i << " is not found in the expression!" << endl;
		}
//		LOG << "A" << i << " is at " << o[i] << endl;
	}	
}


