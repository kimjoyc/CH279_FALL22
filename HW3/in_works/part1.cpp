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
    void printinfo(){
        printf("This shell info: R( %1.2f, %1.2f, %1.2f), with angular momentum: %d, coefficient: %1.2f\n" elem_num, R0(0), R0(1), R0(2), l, alpha, d);
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


class Atom
{
    private:
    Shell sh1;
    Shell sh2;
    Shell sh3;

    public:
    arma::vec make_atom_vec(Shell sh1,Shell sh2,Shell sh3){
        arma::vec 3_prim_gaus={{sh1,sh2,sh3}};
        return 3_prim_gaus;
    }


};

class Molecule{
    private:
    int num_electron;
    int charge;

    public:
    Molecule(int num_electron_input, int input_charge):
    num_electron(num_electron_input), charge(input_charge);

    arma::mat make_square_mat(num_electron){
        return arma::mat Overlap_matrix(num_electron, num_electron,arma::fill::zeros);
    }

};


void ReadShellparameter(Shell & sh1, Shell & sh2, Shell & sh3, string &fname)
{
  ifstream in(fname, ios::in);
  string line1, line2,line3;
  double x0, y0, z0, alpha, d;
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
  Shell sh1,sh2,sh3;

  string fname=argv[1];
  ReadShellparameter(sh1,sh2,sh3,fname);




  return EXIT_SUCCESS;
}