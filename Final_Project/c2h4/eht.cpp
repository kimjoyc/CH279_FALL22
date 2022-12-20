/*

4. Extend your Huckel code to include the pairwise atomic corrections needed to obtain
high accuracy for hydrocarbons (see the paper by Voityuk posted under reading ma-
terial), and implement the gradient. Get some results and compare to your CNDO/2
code for accuracy of structures and relative energies.

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
using namespace arma;


class Shell
{
    private:
    arma::vec R0;
    double alpha;
    double d;
    arma::vec l0;
    int elem_num;
    double gamma_param;
    double beta;

    //constructor
    public:
    Shell(int input_elem_num,double x0_input, double y0_input, double z0_input, double alpha_input, double d_input,double l0_input,double l1_input,double l2_input, double gamma_, double beta_):
    elem_num(input_elem_num),alpha(alpha_input), d(d_input), gamma_param(gamma_), beta(beta_) {l0={l0_input, l1_input, l2_input};R0={x0_input, y0_input, z0_input};} 
    Shell(): alpha(0.5) {l0.zeros(3);R0.zeros(3);}
    ~Shell(){}

    void Reset(int input_elem_num,double x0_input, double y0_input, double z0_input, double alpha_input, double d_input, double l0_input,double l1_input,double l2_input,double gamma_,double beta_)
    {
        elem_num = input_elem_num;R0(0)=x0_input; R0(1)=y0_input; R0(2)=z0_input; alpha=alpha_input; d=d_input; l0(0)=l0_input; l0(1)=l1_input;l0(2)=l2_input; gamma_param=gamma_; beta=beta_;
    }


    arma::vec get_l(){ return l0;}
    double get_alpha(){ return alpha;}
    double get_d(){return d;}
    int get_elem_num(){return elem_num;}
    double get_gamma(){return gamma_param;}
    double get_beta(){return beta;}
    arma::vec get_R0(){ return R0;}
};

void ReadShellparameter(Shell& sh1, Shell& sh2, Shell& sh3, string &fname)
{
  ifstream in(fname, ios::in);
  string line1, line2,line3;
  double x0, y0, z0, alpha, d;
  int elem_num;
  double l0,l1,l2;
  double gamma,beta;
  getline(in, line1);
  getline(in, line2);
  getline(in, line3);

  istringstream iss1(line1);
  if (!(iss1 >> elem_num >> x0 >> y0 >> z0 >> alpha >> d >> l0>>l1>>l2>>gamma>>beta ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh1.Reset(elem_num,x0, y0, z0, alpha,d, l0,l1,l2,gamma,beta);
  istringstream iss2(line2);
  if (!(iss2 >> elem_num >> x0 >> y0 >> z0 >> alpha >> d >> l0>>l1>>l2>>gamma>>beta ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh2.Reset(elem_num,x0, y0, z0, alpha,d, l0,l1,l2,gamma,beta);
  
  istringstream iss3(line3);
  if (!(iss3 >> elem_num >> x0 >> y0 >> z0 >> alpha >> d >> l0>>l1>>l2>>gamma>>beta ))
  {
    throw invalid_argument("There is some problem with format.");
  }
  sh3.Reset(elem_num,x0, y0, z0, alpha,d, l0,l1,l2,gamma,beta);
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

double Eval_Ov(Shell &sh1, Shell & sh2)
{

  arma::vec la = sh1.get_l();
  arma::vec lb = sh2.get_l();

  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();

  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();

  double overlap_tot=Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la(0), lb(0)) 
  * Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la(1), lb(1)) 
  * Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la(2),lb(2));
  
  return overlap_tot;
}


double normalize_func(Shell &sh1)
{
  double same_orb=Eval_Ov(sh1,sh1);
  same_orb=1/sqrt(same_orb);
  return same_orb;
}

double sum_func_opt(Shell* arr1,Shell* arr2)
{
    double ov_tot=0;
    for(int i=0;i<=2;i++)
    {
        for(int k=0;k<=2;k++)
        {
            double orb_1=normalize_func(arr1[i]);

            double orb_2=normalize_func(arr2[k]);

            double ov_elem=Eval_Ov(arr1[i],arr2[k]);

            ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1*orb_2*ov_elem;

   
        }
    }
    
    return ov_tot;    

}

void create_ov_mat(arma::mat &overlap, Shell** arr1)
{
    double row_number=size(overlap)[0];
    double col_number= size(overlap)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {

          overlap(i,j)=sum_func_opt(arr1[i],arr1[j]);

          // if(i==j)
          // {
          //   overlap(i,j)=0;
          // }

          // else
          // {
          //   overlap(i,j)=sum_func_opt(arr1[i],arr1[j]);

          // }

        }
    }

}



void create_ov_mat_offdiag(arma::mat &overlap, Shell** arr1)
{
    double row_number=size(overlap)[0];
    double col_number= size(overlap)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {

          // overlap(i,j)=sum_func_opt(arr1[i],arr1[j]);

          if(i==j)
          {
            overlap(i,j)=0;
          }

          else
          {
            overlap(i,j)=sum_func_opt(arr1[i],arr1[j]);

          }

        }
    }

}

void create_u_mu(arma::mat &overlap,Shell** arr1)
{
    double row_number=size(overlap)[0];
    double col_number= size(overlap)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {


          if(i!=j)
          {
            if(arr1[0][i].get_elem_num()==1)
            {

              overlap(i,j)=-13.605;


            }

            if(arr1[0][i].get_elem_num()==4)
            {
              if(arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0)
              {
                overlap(i,j)=-21.559;
              }

              else
              {
                overlap(i,j)=-13.507;
              }

            }
          }

          else
          {

            overlap(i,j)=0;

          }

        }
    }

}

arma::vec get_g_params(Shell &sh1,Shell&sh2)
{
    double alpha_ab;
    double gamma_ab;
    double omega_ab;
    double r_ab;
    arma::vec h_h_param;
    arma::vec c_c_param;
    arma::vec c_h_param;


    if(sh1.get_elem_num()==sh2.get_elem_num())
    {
        if(sh1.get_elem_num()==1)
        {
            alpha_ab=2.823;
            gamma_ab=12.512;
            omega_ab=-0.0791;
            r_ab=2.279;

            h_h_param={alpha_ab,gamma_ab,omega_ab,r_ab};
            return h_h_param;
        }

        if(sh1.get_elem_num()==4)
        {
            alpha_ab=3.401;
            gamma_ab=658.659;
            omega_ab=0.0312;
            r_ab=3.044;

            c_c_param={alpha_ab,gamma_ab,omega_ab,r_ab};
            return c_c_param;
        }
    }
    
    if(sh1.get_elem_num()!=sh2.get_elem_num())
    {
        if((sh1.get_elem_num()==1&&sh2.get_elem_num()==4)||(sh1.get_elem_num()==4&&sh2.get_elem_num()==1))
        {
            alpha_ab=2.831;
            gamma_ab=99.370;
            omega_ab=-0.0340;
            r_ab=2.843;

            c_h_param={alpha_ab,gamma_ab,omega_ab,r_ab};
            return c_h_param;
            
        }

    }

  return 0;
}

double Solve_EH(arma::mat &OV_mat, arma::mat &H_mat, arma::mat &C_mat, arma::vec &energy_vec, int num_ele){
  //Calculate X_mat
  mat U;
  vec S_eigenvalue;
  arma::eig_sym(S_eigenvalue, U, OV_mat);
  mat S_invsqrt = arma::inv( arma::diagmat( arma::sqrt(S_eigenvalue)) );
  // mat X_mat = U * S_invsqrt;
  mat X_mat = U * S_invsqrt * U.t();
  X_mat.print("X_mat");

  mat H_new = X_mat.t() * H_mat * X_mat;
  // H_new.print("H_new");
  
  arma::mat V;
  arma::eig_sym(energy_vec, V, H_new);
  C_mat = X_mat * V;

  if(num_ele % 2 == 1){
    printf("Warn:: this job is unrestricted");
  }
  double Energy = 2* arma::accu(energy_vec.subvec(0, num_ele/2 - 1));
  return Energy;
}

double h_mu_nu_calc(arma::mat U_mu,arma::mat delta_mu_nu)
{
  double H_mu_nu;
  return H_mu_nu=dot(U_mu,delta_mu_nu);

}
//organize by bond type and then filter by orbital pairs
double h_mu_nu_calc_correction(double H_mu_nu,double beta_mu_nu, double R_ab, double a_o, double lambda_mu_nu)
{
  double param=-lambda_mu_nu*pow(R_ab,2)/pow(a_o,2);
  return H_mu_nu+=beta_mu_nu*sqrt(R_ab/a_o)*exp(param);
  // return H_mu_nu-=beta_mu_nu*sqrt(R_ab/a_o)*exp(param);

}



void create_hmunu(arma::mat &overlap, arma::mat delta_mu_nu,arma::mat U_mu,Shell** arr1)
{
  double row_number=size(overlap)[0];
  double col_number= size(overlap)[1];

  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i==j)
      {
        if(arr1[0][i].get_elem_num()==1)
        {
          overlap(i,j)=-13.605;
        }

        if(arr1[0][i].get_elem_num()==4)
        {
          if(arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0)
          {
            overlap(i,j)=-21.559;
          }

          else
          {
            overlap(i,j)=-13.507;
          }

        }

      }

      else
      {
        create_ov_mat(delta_mu_nu,arr1);
        create_ov_mat_offdiag(U_mu,arr1);
        double h_mu_nu=h_mu_nu_calc(U_mu,delta_mu_nu);
        double lambda_ss;
        double beta_ss;
        double lambda_sp;
        double beta_sp;

        double lambda_pp_sigma;
        double beta_pp_sigma;

        double lambda_pp_pi;
        double beta_pp_pi;


        arma::vec R_ab;
        arma::vec r_ab_sqrd;
        double r_ab_dist;
        double a0=0.52917;


        //only c-h bonds
        if((arr1[0][i].get_elem_num()==1&&arr1[0][j].get_elem_num()==4)||(arr1[0][j].get_elem_num()==4&&arr1[0][i].get_elem_num()==1))
        {
          //ss pair
          if((arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0)&&(arr1[0][j].get_l()[0]==0&&arr1[0][j].get_l()[1]==0&&arr1[0][j].get_l()[2]==0)||
          (arr1[0][j].get_l()[0]==0&&arr1[0][j].get_l()[1]==0&&arr1[0][j].get_l()[2]==0)&&(arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0))
          {
              
            lambda_ss=0.275;
            beta_ss=-8.574;

            R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
            r_ab_sqrd=pow(R_ab,2);
            r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);

            overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_ss,r_ab_dist,a0,lambda_ss);

          }
            //sp
          else
          {
            lambda_sp=0.218;
            beta_sp=-6.813;

            R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
            r_ab_sqrd=pow(R_ab,2);
            r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);

            overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_sp,r_ab_dist,a0,lambda_sp);

          }

        }
        //h-h
        if((arr1[0][i].get_elem_num()==1&&arr1[0][j].get_elem_num()==1)||(arr1[0][j].get_elem_num()==1&&arr1[0][i].get_elem_num()==1))
        {
          lambda_ss=0.280;
          beta_ss=-4.442;

          R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
          r_ab_sqrd=pow(R_ab,2);
          r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);

          overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_ss,r_ab_dist,a0,lambda_ss);
        }

        //c-c
        if(arr1[0][i].get_elem_num()==1&&arr1[0][j].get_elem_num()==1||arr1[0][j].get_elem_num()==1&&arr1[0][i].get_elem_num()==1)
        {
          //ss pair
          if((arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0)&&(arr1[0][j].get_l()[0]==0&&arr1[0][j].get_l()[1]==0&&arr1[0][j].get_l()[2]==0)||
          (arr1[0][j].get_l()[0]==0&&arr1[0][j].get_l()[1]==0&&arr1[0][j].get_l()[2]==0)&&(arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0))
          {
            lambda_ss=0.086;
            beta_ss=-5.969;

            R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
            r_ab_sqrd=pow(R_ab,2);
            r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);
            overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_ss,r_ab_dist,a0,lambda_ss);

          }
          //sp
          if((arr1[0][i].get_l()[0]==1&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0)&&(arr1[0][j].get_l()[0]==0&&arr1[0][j].get_l()[1]==0&&arr1[0][j].get_l()[2]==0)||
          (arr1[0][j].get_l()[0]==1&&arr1[0][j].get_l()[1]==0&&arr1[0][j].get_l()[2]==0)&&(arr1[0][i].get_l()[0]==0&&arr1[0][i].get_l()[1]==0&&arr1[0][i].get_l()[2]==0))
          {
            lambda_sp=0.180;
            beta_sp=-6.160;

            R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
            r_ab_sqrd=pow(R_ab,2);
            r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);

            overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_sp,r_ab_dist,a0,lambda_sp);

          }

          else
          {
            if(arr1[0][i].get_R0()[0]==arr1[0][j].get_R0()[0])
            {
              lambda_pp_sigma=0.186;
              beta_pp_sigma=-8.420;

              R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
              r_ab_sqrd=pow(R_ab,2);
              r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);

              overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_pp_sigma,r_ab_dist,a0,lambda_pp_sigma);

            }
            else
            {
              lambda_pp_pi=0.282;
              beta_pp_pi=-7.403;

              R_ab=arr1[0][i].get_R0()-arr1[0][j].get_R0();
              r_ab_sqrd=pow(R_ab,2);
              r_ab_dist=sqrt(r_ab_sqrd[0]+r_ab_sqrd[1]+r_ab_sqrd[2]);
              overlap(i,j)=h_mu_nu_calc_correction(h_mu_nu,beta_pp_pi,r_ab_dist,a0,lambda_pp_pi);

            }
          }
        }
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


  //ov 
  arma::mat overlap(2,2,arma::fill::zeros);
  Shell* arr1[4]={shell_arr,shell_arr_2,shell_arr_2,shell_arr};
  create_ov_mat(overlap,arr1);
  overlap.print();
  cout << "\n\n";


  arma::mat u_mu(2,2,arma::fill::zeros);
  create_u_mu(u_mu,arr1);
  u_mu.print();
  cout << "\n\n";


  arma::mat overlap_off(2,2,arma::fill::zeros);
  create_ov_mat_offdiag(overlap_off,arr1);
  overlap_off.print();
  cout << "\n\n";

  
  double hmu=h_mu_nu_calc(u_mu,overlap_off);
  cout<<h_mu_nu_calc(u_mu,overlap_off);
  cout << "\n\n";

  // double h_mu_nu_calc_correction(double H_mu_nu,double beta_mu_nu, double R_ab, double a_o, double lambda_mu_nu)
  arma::vec rab=sh1.get_R0()-sh4.get_R0();
  arma::vec rab_=pow(rab,2);
  double rab_dist=sqrt(rab_[0]+rab_[1]+rab_[2]);

  cout << h_mu_nu_calc_correction(hmu,-4.442,rab_dist,0.52917,0.280);
  cout << "\n\n";

  // hmunu
  // void create_hmunu(arma::mat &overlap, arma::mat delta_mu_nu,arma::mat U_mu,Shell** arr1)
  arma::mat hmat(2,2,arma::fill::zeros);
  create_hmunu(hmat,overlap,overlap_off,arr1);
  hmat.print();
  cout << "\n\n";
  return EXIT_SUCCESS;
}