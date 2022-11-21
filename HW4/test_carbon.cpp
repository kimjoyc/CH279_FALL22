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
    Shell(): elem_num(0) {R0.zeros(3);}
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
  sh1.Reset(elem_num,x0, y0, z0, alpha,d, l);
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

arma::mat Eval_Ov(Shell &sh1, Shell & sh2)
{
  int dim1 = sh1.dim_func(), dim2 = sh2.dim_func();
  int la = sh1.get_l(), lb = sh2.get_l();
  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();

  arma::mat overlap_tot(sh1.dim_func(), sh2.dim_func(),arma::fill::zeros);

  int a_index = 0;
  for(int la_x = la; la_x >= 0; la_x--)
    for(int la_y = la - la_x; la_y >= 0; la_y--){
      int b_index = 0;
      for(int lb_x = lb; lb_x >= 0 ; lb_x--)
        for(int lb_y = lb- lb_x; lb_y >= 0; lb_y--) {
          overlap_tot(a_index,b_index)=Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la_x, lb_x) *
              Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la_y, lb_y) * 
              Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la - la_x - la_y, lb -lb_x - lb_y);

          b_index += 1;
        }
      a_index += 1;
    }

    return overlap_tot;

}


arma::mat normalize_func(Shell &sh1){
  arma::mat same_orb=Eval_Ov(sh1,sh1);
  if(sh1.get_l()==1)
  {
    same_orb=1/sqrt(same_orb.diag());
  }

  else 
  {
    same_orb=1/sqrt(same_orb);
  }
  return same_orb;
}

//fix not compatible w non s-orbitals
arma::mat sum_func_opt(Shell* arr1,Shell* arr2)
{
    arma::mat ov_tot(arr1[0].dim_func(), arr2[0].dim_func(),arma::fill::zeros);

    for(int i=0;i<=2;i++)
    {
        for(int k=0;k<=2;k++)
        {
            //pp
            if(arr1[i].get_l()==1&&arr2[k].get_l()==1)
            {
                arma::mat orb_1=normalize_func(arr1[i]);
                arma::mat orb_2=normalize_func(arr2[k]);
                arma::mat ov_elem=Eval_Ov(arr1[i],arr2[k]);
                ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1(0)*orb_2(0)*ov_elem;
            }
            //sp
            if((arr1[i].get_l()==0&&arr2[k].get_l()==1)||(arr1[i].get_l()==1&&arr2[k].get_l()==0))
            {
                arma::mat orb_1=normalize_func(arr1[i]);
                arma::mat orb_2=normalize_func(arr2[k]);
                arma::mat ov_elem=Eval_Ov(arr1[i],arr2[k]);
                ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1(0)*orb_2(0)*ov_elem;
            }
            //ss
            if(arr1[i].get_l()==0&&arr2[k].get_l()==0)
            {
                arma::mat orb_1=normalize_func(arr1[i]);
                arma::mat orb_2=normalize_func(arr2[k]);
                arma::mat ov_elem=Eval_Ov(arr1[i],arr2[k]);
                ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1*orb_2*ov_elem;
            }

        }
    }

    return ov_tot;    

}

double gamma_func(Shell* arr1,Shell* arr2)
{
  double tot=0;
  for(int i=0;i<=2;i++)
  {
    for(int k=0;k<=2;k++)
    {
      for(int l=0;l<=2;l++)
      {
        for(int m=0;m<=2;m++)
        {
          arma::mat orb_1=normalize_func(arr1[i]);
          arma::mat orb_2=normalize_func(arr1[k]);
          arma::mat orb_3=normalize_func(arr2[l]);
          arma::mat orb_4=normalize_func(arr2[m]);

          arma::vec Ra=arr1[i].get_R0();

          
          arma::vec Rb=arr2[i].get_R0();
          arma::vec r_ab=Ra-Rb;
          arma::vec r_ab_=pow(r_ab,2);
          double R_ab_dist=sqrt(r_ab_(0)+r_ab_(1)+r_ab_(2));

        
          double omeg_A=1/(arr1[i].get_alpha()+arr1[k].get_alpha());
          double omeg_B=1/(arr2[l].get_alpha()+arr2[m].get_alpha());


          double Ua=powf(M_PI*omeg_A,1.5);
          double Ub=powf(M_PI*omeg_B,1.5);
          double V_squared=1/(omeg_A+omeg_B);

          arma::vec T=V_squared*pow(Ra-Rb,2);

          double T_=sqrt(T(0)+T(1)+T(2));
          double U=Ua*Ub;
          double zer_AB=U*erf(T_)*1/R_ab_dist;


          double d_prime_1=arr1[i].get_d()*orb_1(0);
          double d_prime_2=arr1[k].get_d()*orb_2(0);
          double d_prime_3=arr2[l].get_d()*orb_3(0);
          double d_prime_4=arr2[m].get_d()*orb_4(0);

          double zer_AA=U*sqrt(2*V_squared)*sqrt(2/M_PI);

          if(R_ab_dist==0)
          {
            tot+=d_prime_1*d_prime_2*d_prime_3*d_prime_4*zer_AA;
          }
          else
          {
            tot+=d_prime_1*d_prime_2*d_prime_3*d_prime_4*zer_AB;
          }


        }
      }
    }
  }
  return tot*27.2114;
}


double h_core_diag(double semi_emp, double Z_a, double Z_b, double gamma_AA, double gamma_AB){

  double h_core = -semi_emp-(Z_a-0.5)*gamma_AA - (Z_b*gamma_AB);
  return h_core;
}

double h_core_diag_off(double beta_a, double beta_b,double S_mu_nu)
{
  double h_core = 0.5*(beta_a+beta_b)*S_mu_nu;
  return h_core;
}


arma::mat scf(arma::mat h_core){
  arma::mat F_a=h_core;
  arma::mat F_b=h_core;
  arma::vec epsilon_a;
  arma::mat C_a;
  eig_sym(epsilon_a,C_a,F_a);


  return C_a.col(0)*C_a.col(0).t();
}


double g_mat_off(double P_mat,double gamma_AB){
  return -P_mat*gamma_AB;

}

double g_mat(double P_aa_tot, double P_mu_mu, double P_bb_tot,double gamma_AA,double gamma_AB)
{
  return (P_aa_tot-P_mu_mu)*gamma_AA+P_bb_tot*gamma_AB;

}


double scf_real_v_nuc(double Z_a, double Z_b, double n_atom, Shell* arr1, Shell* arr2)
{
  double V_nuc=0;
  for (int i = 0; i < n_atom; i++)
  {
    for (int j = 0; j < i; j++) 
    {
      arma::vec Ra=arr1[0].get_R0();
      arma::vec Rb=arr2[0].get_R0();
      arma::vec r_ab=Ra-Rb;    
      arma::vec r_ab_=pow(r_ab,2);
      double R_ab_dist=sqrt(r_ab_(0)+r_ab_(1)+r_ab_(2));
      V_nuc+=(Z_a*Z_b)/R_ab_dist;

    }
  }
  return V_nuc*27.2114;

}


void create_gamma_mat(arma::mat &gamma, Shell** arr1)
{
    double row_number=size(gamma)[0];
    double col_number= size(gamma)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
            gamma(i,j)=gamma_func(arr1[i],arr1[j]);
        }
    }
}

void create_ov_mat(arma::mat overlap, Shell** arr1)
{
    double row_number=size(overlap)[0];
    double col_number= size(overlap)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
            if(arr1[i][0].get_l()==1&&arr1[j][0].get_l()==1)
            {
                overlap(i,j)=sum_func_opt(arr1[i],arr1[j]).diag()[0];
                overlap(i,j)=sum_func_opt(arr1[i],arr1[j]).diag()[1];
                overlap(i,j)=sum_func_opt(arr1[i],arr1[j]).diag()[2];

            }


        }
    }
}

int main(int argc, char* argv[])
{
//H1 
  Shell sh1,sh2,sh3;
  string fname_1=argv[1];
  ReadShellparameter(sh1,sh2,sh3,fname_1);
  Shell shell_arr[3]={sh1,sh2,sh3};
//H3 
  Shell sh4,sh5,sh6;
  string fname_2=argv[2];
  ReadShellparameter(sh4,sh5,sh6,fname_2);
  Shell shell_arr_2[3]={sh4,sh5,sh6};

//C1-2s
  Shell sh7,sh8,sh9;
  string fname_3=argv[3];
  ReadShellparameter(sh7,sh8,sh9,fname_3);
  Shell shell_arr_3[3]={sh7,sh8,sh9};
  
//C1-2p
  Shell sh10,sh11,sh12;
  string fname_4=argv[4];
  ReadShellparameter(sh10,sh11,sh12,fname_4);
  Shell shell_arr_4[3]={sh10,sh11,sh12};

//C2-2s
  Shell sh13,sh14,sh15;
  string fname_5=argv[5];
  ReadShellparameter(sh13,sh14,sh15,fname_5);
  Shell shell_arr_5[3]={sh13,sh14,sh15};
//C2-2p
  Shell sh16,sh17,sh18;
  string fname_6=argv[6];
  ReadShellparameter(sh16,sh17,sh18,fname_6);
  Shell shell_arr_6[3]={sh16,sh17,sh18};

  //gamma
  arma::mat gamma(4, 4,arma::fill::zeros);
  Shell* arr1[4]={shell_arr,shell_arr_3,shell_arr_5,shell_arr_2};
  create_gamma_mat(gamma,arr1);
//   gamma.print();

  //ov
  arma::mat overlap(10,10,arma::fill::zeros);
  Shell* arr2[10]={shell_arr,shell_arr_3,shell_arr_4,shell_arr_4,shell_arr_4,shell_arr_5,shell_arr_6,shell_arr_6,shell_arr_6,shell_arr_2};
  create_ov_mat(overlap,arr2);
  overlap.print();


  
  return EXIT_SUCCESS;
}