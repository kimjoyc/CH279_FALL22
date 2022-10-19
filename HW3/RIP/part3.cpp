/*
Your program should perform the following tasks:
• Read in the coordinates, in the format: E X Y Z for each atom, where E is the element
(handle at least H and C).
• Evaluate the number of basis functions, N from the molecular formula, Ca Hb , where
the relation is N =4a +b. Your matrices, such as S, H, X, etc, will be N ×N , so this will
enable you to define them.


Evaluate the number of electrons 2n =4a +b. Throw an error if the number of electron
pairs n =2a +b/2 is not an integer. Knowing n is necessary to evaluate the energy later.
• Build a list of the basis functions, which are contracted gaussians, as described in de-
tail below. For each basis function, ωμ(r) (for μ=1 ···N ) there will be (i) a center, R, (ii)
3 quantum numbers, (l , m, n), and (iii) information about 3 primitive functions: 3 ex-
ponents, αk (given below), 3 corresponding contraction coefficients, dk (given below),
and 3 normalization constants, N l mn
k (to be worked out by your code in q. 2 below).
See Eqs. 1.1 and 1.2 below for definitions.
As a reminder, the atomic orbital (AO) basis functio



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
    int elem_num;
    int N;
    int electro;

    //constructor
    public:
    Shell(int elem_num_input, int elec_input, double x0_input, double y0_input, double z0_input, double alpha_input, double d_input, int l_input):
    electro(elec_input),elem_num(elem_num_input),alpha(alpha_input), d(d_input), l(l_input) {R0={x0_input, y0_input, z0_input};}
    Shell():alpha(0.5), l(0) {R0.zeros(3);}
    ~Shell(){}


    void Reset(double x0_input, double y0_input, double z0_input){
        R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input; 
    }

    void Reset2(double alpha_input, double d_input, int l_input){
        alpha=alpha_input; d=d_input; l=l_input;

    }

    void Reset3(int input_N, int elec_input){
        N=input_N; electro=elec_input;

    }

    void printinfo(){
        printf("This shell info: R( %1.2f, %1.2f, %1.2f), with angular momentum: %d, coefficient: %1.2f\n",R0(0), R0(1), R0(2));
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


void ReadShellparameter(Shell & sh1, Shell & sh2, Shell & sh3, string &fname)
{
  ifstream in(fname, ios::in);
  string line1, line2,line3;
  int N;
  int electro;
  int elem_num; 
  double x0, y0, z0;
  getline(in, line1);
  getline(in, line2);
  getline(in, line3);

  istringstream iss1(line1);
  if (!(iss1 >> N>>electro ))
  {
    throw invalid_argument("Thermat.");
  }
  sh1.Reset3(N,electro);
  in.close();


  istringstream iss2(line2);
  if (!(iss2 >> elem_num>> x0 >> y0 >> z0 ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh2.Reset(x0, y0, z0);
  in.close();

  istringstream iss3(line3);
  if (!(iss3 >> elem_num>> x0 >> y0 >> z0 ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh3.Reset(x0, y0, z0);
  in.close();
}

void ReadShell_exp_contraction(Shell & sh1, Shell & sh2, Shell &sh3, string &fname)
{
  ifstream in(fname, ios::in);
  string line1, line2, line3;
  double alpha, d;
  int l;
  getline(in, line1);
  getline(in, line2);
  getline(in, line3);

  istringstream iss1(line1);
  if (!(iss1 >>  alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh1.Reset2(alpha,d, l);
  istringstream iss2(line2);
  if (!(iss2 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh2.Reset2(alpha, d,l);
  in.close();

  istringstream iss3(line3);
  if (!(iss3 >> alpha >> d >> l ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh3.Reset2(alpha, d,l);
  in.close();
}



int main(int argc, char* argv[])
{
  //E X Y Z for each atom, where E is the element (handle at least H and C).
  Shell sh1,sh2,sh3;
  string fname(argv[1]);
  ReadShellparameter(sh1, sh2, sh3,fname);

  cout << sh2.get_R0();

  // //N =4a +b. Your matrices, such as S, H, X, etc, will be N ×N , so this will enable you to define them.
  // string fname_2(argv[2]);
  // ReadShell_exp_contraction(sh1,sh2,sh3,fname_2);




  return EXIT_SUCCESS;
}
