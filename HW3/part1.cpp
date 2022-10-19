/*

*/


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

    //constructor
    public:
    Shell(double x0_input, double y0_input, double z0_input, double alpha_input, double d_input, int l_input):
    alpha(alpha_input), d(d_input), l(l_input) {R0={x0_input, y0_input, z0_input};}
    Shell():alpha(0.5), l(0) {R0.zeros(3);}
    ~Shell(){}


    void Reset(double x0_input, double y0_input, double z0_input, double alpha_input, double d_input, int l_input){
        R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input; alpha=alpha_input; d=d_input; l=l_input;

    }
    void printinfo(){
        printf("This shell info: R( %1.2f, %1.2f, %1.2f), with angular momentum: %d, coefficient: %1.2f\n",R0(0), R0(1), R0(2), l, alpha, d);
    }


    // give the numbers of functions in this shell
    int dim_func()
    {
        return (l+1)*(l+2)/2;

    }
    int get_l(){ return l;}
    double get_alpha(){ return alpha;}
    double get_d(){return d;}
    arma::vec get_R0(){ return R0;}
};


void ReadShellparameter(Shell & sh1, Shell & sh2, string &fname)
{
  ifstream in(fname, ios::in);
  string line1, line2;
  double x0, y0, z0, alpha, d;
  int l;
  getline(in, line1);
  getline(in, line2);
  istringstream iss1(line1);
  if (!(iss1 >> x0 >> y0 >> z0 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh1.Reset(x0, y0, z0, alpha,d, l);
  istringstream iss2(line2);
  if (!(iss2 >> x0 >> y0 >> z0 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh2.Reset(x0, y0, z0, alpha, d,l);
  in.close();
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

void Eval_Ov(arma::mat &Overlap, Shell& sh1, Shell& sh2){
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



int main(int argc, char* argv[])
{
  Shell sh1,sh2;
  int dim1_to = sh1.dim_func(), dim2_to = sh2.dim_func();
  arma::mat Overlap_matrix_tot(dim1_to, dim2_to,arma::fill::zeros);

  for(int i = 1; i < argc-1; i++)
  {
    Shell sh1,sh2;
    string fname(argv[i]);
    int dim1 = sh1.dim_func(), dim2= sh2.dim_func();

    arma::mat Overlap_matrix(dim1,dim2,arma::fill::zeros);
    ReadShellparameter(sh1, sh2, fname);
    Eval_Ov(Overlap_matrix,sh1,sh2);
    Overlap_matrix_tot+=Overlap_matrix*sh1.get_d()*sh2.get_d();
  }

  Overlap_matrix_tot.print();


  return EXIT_SUCCESS;
}
