#include <iostream>
#include <armadillo>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <cassert>
#include "util.h"
#include <stdlib.h>
#include <stdexcept>
#include <armadillo>
#include <stdio.h>

using namespace std;


class Shell
{
    private:
    arma::vec R0;
    double alpha;
    double d;
    int l;
    int elem_num;

    //constructor
    public:
    Shell(int input_elem_num,double x0_input, double y0_input, double z0_input, double alpha_input, double d_input, int l_input):
    elem_num(input_elem_num),alpha(alpha_input), d(d_input), l(l_input) {R0={x0_input, y0_input, z0_input};}
    Shell():alpha(0.5), l(0) {R0.zeros(3);}
    ~Shell(){}


    void Reset(int input_elem_num,double x0_input, double y0_input, double z0_input, double alpha_input, double d_input, int l_input){
        elem_num = input_elem_num;R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input; alpha=alpha_input; d=d_input; l=l_input;

    }


    // give the numbers of functions in this shell
    int dim_func()
    {
        return (l+1)*(l+2)/2;

    }
    int get_l(){ return l;}
    double get_alpha(){ return alpha;}
    double get_d(){return d;}
    int get_elem_num(){return elem_num;}
    arma::vec get_R0(){ return R0;}
};





void ReadShellparameter(Shell& sh1, Shell& sh2, Shell& sh3, string &fname)
{
  ifstream in(fname, ios::in);
  string line1, line2,line3;
  double x0, y0, z0, alpha;
  float d;
  int elem_num;
  int l;
  getline(in, line1);
  getline(in, line2);
  getline(in, line3);

  istringstream iss1(line1);
  if (!(iss1 >> elem_num >> x0 >> y0 >> z0 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  istringstream iss2(line2);
  if (!(iss2 >> elem_num >> x0 >> y0 >> z0 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh2.Reset(elem_num,x0, y0, z0, alpha,d, l);
  
  istringstream iss3(line3);
  if (!(iss3 >> elem_num >> x0 >> y0 >> z0 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh3.Reset(elem_num,x0, y0, z0, alpha,d, l);
}

double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb)
{
  double prefactor = exp( -alphaa*alphab*(xa-xb)*(xa-xb)/(alphaa+ alphab))* sqrt(M_PI / (alphaa+ alphab)) ;
  double xP = (alphaa* xa + alphab * xb)/ (alphaa+ alphab);

  double result = 0.0;
  for(int i_index = 0; i_index <= la; i_index++)
    for(int j_index = 0; j_index <= lb; j_index++){
      if((i_index + j_index) % 2 == 1)
        continue;
      double C_part = Combination(la, i_index) * Combination(lb, j_index);
      double DF_part = DoubleFactorial(i_index + j_index - 1);
      double numerator = pow(xP-xa, la -i_index) * pow(xP-xb, lb - j_index);
      double denominator = pow(2*(alphaa+ alphab), double(i_index + j_index ) / 2.0);
      double temp = C_part * DF_part * numerator / denominator;
      result += temp;
    }

  result *= prefactor;
  return result;
}

void Eval_Ov(arma::mat & Overlap, Shell &sh1, Shell & sh2){
  int dim1 = sh1.dim_func(), dim2 = sh2.dim_func();
  assert(Overlap.n_rows == dim1);
  assert(Overlap.n_cols == dim2);

  int la = sh1.get_l(), lb = sh2.get_l();
  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();

  int a_index = 0;
  for(int la_x = la; la_x >= 0; la_x--)
    for(int la_y = la - la_x; la_y >= 0; la_y--){
      int b_index = 0;
      for(int lb_x = lb; lb_x >= 0 ; lb_x--)
        for(int lb_y = lb- lb_x; lb_y >= 0; lb_y--) {
          Overlap(a_index, b_index) = Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la_x, lb_x) *
              Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la_y, lb_y) * 
              Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la - la_x - la_y, lb -lb_x - lb_y);

          b_index += 1;
        }
      a_index += 1;
    }

}

void normalize_func(arma::mat &same_orb,Shell &sh1){
  Eval_Ov(same_orb, sh1,sh1);
  same_orb=1/sqrt(same_orb);
}


// void sum_func(arma::mat & ov_tot, Shell &sh1,Shell &sh2, Shell &sh3, Shell &sh4, Shell &sh5, Shell &sh6){
//   for(int i=0;i<=3;i++)
//   {
//     for(int k=0;k<=3;k++)
//     {
//       arma::mat orb_i(sh1.dim_func(), sh1.dim_func(),arma::fill::zeros);
//       normalize_func(orb_1,sh1);
//       arma::mat orb_k(sh2.dim_func(), sh2.dim_func(),arma::fill::zeros);
//       normalize_func(orb_2,sh2);
//       arma::mat orb_k(sh3.dim_func(), sh3.dim_func(),arma::fill::zeros);
//       normalize_func(orb_3,sh3);

//       arma::mat orb_k(sh4.dim_func(), sh4.dim_func(),arma::fill::zeros);
//       normalize_func(orb_4,sh4);
//       arma::mat orb_k(sh5.dim_func(), sh5.dim_func(),arma::fill::zeros);
//       normalize_func(orb_5,sh5);
//       arma::mat orb_k(sh6.dim_func(), sh6.dim_func(),arma::fill::zeros);
//       normalize_func(orb_6,sh6);


//       arma::mat Overlap_matrix(sh1.dim_func(), sh2.dim_func(),arma::fill::zeros);
//       Eval_Ov(Overlap_matrix,sh1,sh2);

//       ov_tot=sh1.get_d()*sh2.get_d()*orb_i*orb_k*Overlap_matrix;

//     }

//   }
    
// }

//arr of ptrs

void sum_func(arma::mat & ov_tot, Shell &sh1,Shell &sh2, Shell &sh3, Shell &sh4, Shell &sh5, Shell &sh6){
  
  //atom1
  arma::mat orb_1(sh1.dim_func(), sh1.dim_func(),arma::fill::zeros);
  normalize_func(orb_1,sh1);
  arma::mat orb_2(sh2.dim_func(), sh2.dim_func(),arma::fill::zeros);
  normalize_func(orb_2,sh2);
  arma::mat orb_3(sh3.dim_func(), sh3.dim_func(),arma::fill::zeros);
  normalize_func(orb_3,sh3);
  //atom2
  arma::mat orb_4(sh4.dim_func(), sh4.dim_func(),arma::fill::zeros);
  normalize_func(orb_4,sh4);
  arma::mat orb_5(sh5.dim_func(), sh5.dim_func(),arma::fill::zeros);
  normalize_func(orb_5,sh5);
  arma::mat orb_6(sh6.dim_func(), sh6.dim_func(),arma::fill::zeros);
  normalize_func(orb_6,sh6);


  //atom1 first gauss with atom2 3 gauss
  arma::mat Overlap_matrix_1_4(sh1.dim_func(), sh4.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_1_4,sh1,sh4);
  arma::mat Overlap_matrix_1_5(sh1.dim_func(), sh5.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_1_5,sh1,sh5);
  arma::mat Overlap_matrix_1_6(sh1.dim_func(), sh6.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_1_6,sh1,sh6);

  ov_tot+=sh1.get_d()*sh4.get_d()*orb_1*orb_4*Overlap_matrix_1_4;
  ov_tot+=sh1.get_d()*sh5.get_d()*orb_1*orb_5*Overlap_matrix_1_5;
  ov_tot+=sh1.get_d()*sh6.get_d()*orb_1*orb_6*Overlap_matrix_1_6;

  //atom1 2nd gauss w/atom2 3 gauss
  arma::mat Overlap_matrix_2_4(sh2.dim_func(), sh4.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_2_4,sh2,sh4);
  arma::mat Overlap_matrix_2_5(sh2.dim_func(), sh5.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_2_5,sh2,sh5);
  arma::mat Overlap_matrix_2_6(sh2.dim_func(), sh6.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_2_6,sh2,sh6);

  ov_tot+=sh2.get_d()*sh4.get_d()*orb_2*orb_4*Overlap_matrix_2_4;
  ov_tot+=sh2.get_d()*sh5.get_d()*orb_2*orb_5*Overlap_matrix_2_5;
  ov_tot+=sh2.get_d()*sh6.get_d()*orb_2*orb_6*Overlap_matrix_2_6;

  //atom1 3rd gauss w/atom2 3 gauss
  arma::mat Overlap_matrix_3_4(sh3.dim_func(), sh4.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_3_4,sh3,sh4);
  arma::mat Overlap_matrix_3_5(sh3.dim_func(), sh5.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_3_5,sh3,sh5);
  arma::mat Overlap_matrix_3_6(sh3.dim_func(), sh6.dim_func(),arma::fill::zeros);
  Eval_Ov(Overlap_matrix_3_6,sh3,sh6);

  ov_tot+=sh3.get_d()*sh4.get_d()*orb_3*orb_4*Overlap_matrix_3_4;
  ov_tot+=sh3.get_d()*sh5.get_d()*orb_3*orb_5*Overlap_matrix_3_5;
  ov_tot+=sh3.get_d()*sh6.get_d()*orb_3*orb_6*Overlap_matrix_3_6;

}


//matrix manipulation fn

//1. make the orthogonalization transformation: X =S−1/2
void orth_trans(arma::mat &Overlap){
    
    arma::vec S_eigval;
    arma::mat S_eigvec;
    eig_sym(S_eigval, S_eigvec, Overlap);  
    Overlap= S_eigvec*sqrt(inv(diagmat(S_eigval)))*S_eigvec.t();
}

// X = S^1/2
void orth_trans_(arma::mat &Overlap){
    
    arma::vec S_eigval;
    arma::mat S_eigvec;
    eig_sym(S_eigval, S_eigvec, Overlap);  
    Overlap= S_eigvec*sqrt(diagmat(S_eigval))*S_eigvec.t();
}

//ham matrix 
void h_mat_get(arma::mat &h2,arma::mat &h){
    h2=1.75*0.5*arma::sum(h2.diag(),0);
    h=h2[0]*h;
}


//. form the hamiltonian in the orthogonalized basis: H=XT HX
void ham_ortho(arma::mat &Hmat,arma::mat &S_inv){
    Hmat= S_inv.t()*Hmat*S_inv;
}

//3. diagonalize: HV =Vε
void h_diag(arma::mat &h, arma::mat &s){
    // S_inv.t()*Hmat*S_inv
    arma::vec epsilon;
    arma::mat v_prim;
    eig_sym(epsilon, v_prim, h);  

    // s = s*v_prim;
}


//4. form the MO coefficients: C =XV 
//basis vectors form an orthonormal set, the orbital overlap matrix will be the identity matrix. T
void get_v(arma::mat &H1, arma::mat &H2){

  orth_trans(H1);
  orth_trans_(H2);
  H1=H1*H2;
}

//calc total energy
/*
You evaluate this by
subtracting twice the energy of the H atom (−27.2 eV ) from your computed H2 energy.
This result is also your next debugging case. The result should be around −4 eV (the
minus sign indicates binding).
*/

void total_energy(arma::mat &h){
    
    arma::vec epsilon;
    arma::mat v_prim;
    eig_sym(epsilon, v_prim, h);    
    h= arma::sum(epsilon,0)*2;
}

void total_energy_sep(arma::mat &h){
    
    arma::vec epsilon;
    arma::mat v_prim;
    eig_sym(epsilon, v_prim, h);    
    h=epsilon*2;
}

void diff_calc(arma::mat &h2, arma::mat &h1){
    h2= h2[0];
    total_energy(h1);
    h2=h2-h1; 
}





int main(int argc, char* argv[])
{
  Shell sh1,sh2,sh3;
  string fname_1=argv[1];
  ReadShellparameter(sh1,sh2,sh3,fname_1);


  Shell sh4,sh5,sh6;
  string fname_2=argv[2];
  ReadShellparameter(sh4,sh5,sh6,fname_2);

  //h1*h1
  arma::mat ov_tot_1(sh1.dim_func(), sh1.dim_func(),arma::fill::zeros);
  sum_func(ov_tot_1,sh4,sh5,sh6,sh4,sh5,sh6);
  //h2*h1
  arma::mat ov_tot_2(sh1.dim_func(), sh1.dim_func(),arma::fill::zeros);
  sum_func(ov_tot_2,sh1,sh2,sh3,sh4,sh5,sh6);
  arma::mat ov_tot_={{ov_tot_1[0],ov_tot_2[0]},{ov_tot_2[0],ov_tot_1[0]}};

  //h matrix
  arma::mat h_diag= {{-13.6,0},{0,-13.6}};
  h_mat_get(h_diag,ov_tot_);
  arma::mat h_mat = {{-13.6,ov_tot_[1]},{ov_tot_[1],-13.6}};
  h_mat.print();

  // form the hamiltonian in the orthogonalized basis: H=XT HX
  arma::mat ov_tot__={{ov_tot_1[0],ov_tot_2[0]},{ov_tot_2[0],ov_tot_1[0]}};
  //S^-1/2
  orth_trans(ov_tot__);
  ov_tot__.print();
  ham_ortho(ov_tot_,ov_tot__);
  ov_tot_.print();

  // // 3. diagonalize: HV =Vε
  // h_diag(ov_tot_,ov_tot__);
  // ov_tot__.print();



  return EXIT_SUCCESS;
}