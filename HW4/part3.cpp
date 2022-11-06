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


double Eval_Ov(Shell &sh1, Shell & sh2){
  int dim1 = sh1.dim_func(), dim2 = sh2.dim_func();

  int la = sh1.get_l(), lb = sh2.get_l();
  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();


  double overlap_tot=0;

  int a_index = 0;
  for(int la_x = la; la_x >= 0; la_x--)
    for(int la_y = la - la_x; la_y >= 0; la_y--){
      int b_index = 0;
      for(int lb_x = lb; lb_x >= 0 ; lb_x--)
        for(int lb_y = lb- lb_x; lb_y >= 0; lb_y--) {
          overlap_tot+= Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la_x, lb_x) *
              Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la_y, lb_y) * 
              Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la - la_x - la_y, lb -lb_x - lb_y);

          b_index += 1;
        }
      a_index += 1;
    }

    return overlap_tot;

}

double normalize_func(Shell &sh1){
  double same_orb=Eval_Ov(sh1,sh1);
  same_orb=1/sqrt(same_orb);
  return same_orb;
}


double sum_func_opt(Shell* arr1,Shell* arr2){
    double ov_tot=0;
    for(int i=0;i<=2;i++){
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
          double orb_1=normalize_func(arr1[i]);
          double orb_2=normalize_func(arr1[k]);
          double orb_3=normalize_func(arr2[l]);
          double orb_4=normalize_func(arr2[m]);

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


          double d_prime_1=arr1[i].get_d()*orb_1;
          double d_prime_2=arr1[k].get_d()*orb_2;
          double d_prime_3=arr2[l].get_d()*orb_3;
          double d_prime_4=arr2[m].get_d()*orb_4;

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

  // cout << epsilon_a;
  // cout << "\n\n";
  // cout << C_a;
  // cout << "\n\n";
  // cout << C_a.col(0)*C_a.col(0).t();
  // cout << "\n\n";



  return C_a.col(0)*C_a.col(0).t();
}


double g_mat_off(double P_mat,double gamma_AB){
  return -P_mat*gamma_AB;

}

double g_mat(double P_aa_tot, double P_mu_mu, double P_bb_tot,double gamma_AA,double gamma_AB)
{
  return (P_aa_tot-P_mu_mu)*gamma_AA+P_bb_tot*gamma_AB;

}


double scf_real_v_nuc(double Z_a, double Z_b, Shell* arr1, Shell* arr2)
{
  double V_nuc=0;
  for(int i=0;i<=2;i++)
  {
    arma::vec Ra=arr1[i].get_R0();
    // cout << Ra;
    // cout << "\n\n";
    arma::vec Rb=arr2[i].get_R0();
    // cout << Rb;
    // cout << "\n\n";
    arma::vec r_ab=Ra-Rb;    
    // cout << r_ab;
    // cout << "\n\n";
    arma::vec r_ab_=pow(r_ab,2);
    // cout << r_ab_;
    // cout << "\n\n";
    double R_ab_dist=sqrt(r_ab_(0)+r_ab_(1)+r_ab_(2));
    // cout << R_ab_dist;
    // cout << "\n\n";
    V_nuc+=(Z_a*Z_b)/R_ab_dist;
    // cout << V_nuc;
    // cout << "\n\n";
  }

  return V_nuc;

  // double energy_cdno2;
  // double energy_cdno2_new;



  // while(energy_cdno2==energy_cdno2_new)
  // {
  //   arma::mat fock_mat_a = h_mu_nu+g_mu_nu;
  //   arma::mat fock_mat_b = h_mu_nu+g_mu_nu;

  //   arma::vec ep_a;
  //   arma::mat rho_mat_alpha;
  //   eig_sym(ep_a,rho_mat_alpha,fock_mat_a);

  //   arma::vec ep_b;
  //   arma::mat rho_mat_beta;
  //   eig_sym(ep_b,rho_mat_beta,fock_mat_b);

  //   arma::mat C_a=rho_mat_alpha.t()*rho_mat_alpha;
  //   arma::mat C_b=rho_mat_beta.t()*rho_mat_beta;

  //   arma::mat P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
  //   arma::mat P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

  //   arma::mat P_t=P_a_new.col(0)+P_b_new.col(0);
    
  //   double gamma_AA=gamma_func(arr1,arr1);
  //   double gamma_AB=gamma_func(arr2,arr1);

  //   double Pmat_a=P_a_new(0);
  //   double P_aa_tot=sum(P_a_new.col(0));
  //   double P_bb_tot=sum(P_b_new.col(0));
  //   double Pmat_b=P_b_new(0);
    
  //   arma::mat G_a={{g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_a,gamma_AB)},{g_mat_off(Pmat_a,gamma_AB),g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB)}};
  //   arma::mat G_b={{g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_b,gamma_AB)},{g_mat_off(Pmat_b,gamma_AB),g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB)}};

  //   arma::mat F_a=h_mu_nu+G_a; 
  //   arma::mat F_b=h_mu_nu+G_b;


  //   energy_cdno2=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;

  //   //next iter

  //   fock_mat_a =F_a;
  //   fock_mat_b=F_b;

  //   arma::vec ep_a;
  //   arma::mat rho_mat_alpha;
  //   eig_sym(ep_a,rho_mat_alpha,fock_mat_a);

  //   arma::vec ep_b;
  //   arma::mat rho_mat_alpha;
  //   eig_sym(ep_b,rho_mat_beta,fock_mat_b);

  //   C_a=rho_mat_alpha.t()*rho_mat_alpha;
  //   C_b=rho_mat_beta.t()*rho_mat_beta;

  //   P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
  //   P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

  //   P_t=P_a_new.col(0)+P_b_new.col(0);
    
  //   Pmat_a=P_a_new(0);
  //   P_aa_tot=sum(P_a_new.col(0));
  //   P_bb_tot=sum(P_b_new.col(0));
  //   Pmat_b=P_b_new(0);
    
  //   G_a={{g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_a,gamma_AB)},{g_mat_off(Pmat_a,gamma_AB),g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB)}};
  //   G_b={{g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_b,gamma_AB)},{g_mat_off(Pmat_b,gamma_AB),g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB)}};

  //   F_a=h_mu_nu+G_a; 
  //   F_b=h_mu_nu+G_b;

  //   energy_cdno2_new=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;
  // }
}


double energy_cdno2_scf(double V_nuc,double gamma_AA,double gamma_AB,arma::mat h_mu_nu,arma::mat g_mu_nu)
{
  double energy_cdno2;

  arma::mat fock_mat_a = h_mu_nu+g_mu_nu;
  arma::mat fock_mat_b = h_mu_nu+g_mu_nu;
  arma::vec ep_a;
  arma::mat rho_mat_alpha;
  eig_sym(ep_a,rho_mat_alpha,fock_mat_a);
  
  arma::vec ep_b;
  arma::mat rho_mat_beta;
  eig_sym(ep_b,rho_mat_beta,fock_mat_b);

  arma::mat C_a=rho_mat_alpha.t()*rho_mat_alpha;
  arma::mat C_b=rho_mat_beta.t()*rho_mat_beta;

  arma::mat P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
  arma::mat P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

  arma::mat P_t=P_a_new.col(0)+P_b_new.col(0);

  double Pmat_a=P_a_new(0);
  double P_aa_tot=sum(P_a_new.col(0));
  double P_bb_tot=sum(P_b_new.col(0));
  double Pmat_b=P_b_new(0);
    
  arma::mat G_a={{g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_a,gamma_AB)},{g_mat_off(Pmat_a,gamma_AB),g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB)}};
  arma::mat G_b={{g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_b,gamma_AB)},{g_mat_off(Pmat_b,gamma_AB),g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB)}};

  arma::mat F_a=h_mu_nu+G_a; 
  arma::mat F_b=h_mu_nu+G_b;

  energy_cdno2=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;

  return energy_cdno2;


}


// double real_scf(double energy_cdno2,)
// {
//   double energy_cdno2_new;
//   while(energy_cdno2==energy_cdno2_new)
//   {
    

//   }



// }


double scf_real(double V_nuc,double gamma_AA,double gamma_AB,arma::mat h_mu_nu,arma::mat g_mu_nu)
{
  double energy_cdno2;
  double energy_cdno2_new;


  while(energy_cdno2==energy_cdno2_new)
  {
    arma::mat fock_mat_a = h_mu_nu+g_mu_nu;
    arma::mat fock_mat_b = h_mu_nu+g_mu_nu;

    arma::vec ep_a;
    arma::mat rho_mat_alpha;
    eig_sym(ep_a,rho_mat_alpha,fock_mat_a);

    arma::vec ep_b;
    arma::mat rho_mat_beta;
    eig_sym(ep_b,rho_mat_beta,fock_mat_b);

    arma::mat C_a=rho_mat_alpha.t()*rho_mat_alpha;
    arma::mat C_b=rho_mat_beta.t()*rho_mat_beta;

    arma::mat P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
    arma::mat P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

    arma::mat P_t=P_a_new.col(0)+P_b_new.col(0);

    double Pmat_a=P_a_new(0);
    double P_aa_tot=sum(P_a_new.col(0));
    double P_bb_tot=sum(P_b_new.col(0));
    double Pmat_b=P_b_new(0);
    
    arma::mat G_a={{g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_a,gamma_AB)},{g_mat_off(Pmat_a,gamma_AB),g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB)}};
    cout << G_a;
    cout << "\n\n";
    arma::mat G_b={{g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_b,gamma_AB)},{g_mat_off(Pmat_b,gamma_AB),g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB)}};
    cout << G_b;
    cout << "\n\n";

    arma::mat F_a=h_mu_nu+G_a; 
    cout << F_a;
    cout << "\n\n";

    arma::mat F_b=h_mu_nu+G_b;
    cout << F_b;
    cout << "\n\n";

    energy_cdno2=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;
    cout << energy_cdno2;
    cout << "\n\n";

    //next iter
    fock_mat_a=F_a;
    fock_mat_b=F_b;

    eig_sym(ep_a,rho_mat_alpha,fock_mat_a);
    cout << ep_a;
    cout << "\n\n";
    cout<<rho_mat_alpha;
    cout << "\n\n";

    eig_sym(ep_b,rho_mat_beta,fock_mat_b);
    cout << ep_b;
    cout << "\n\n";
    cout<<rho_mat_beta;
    cout << "\n\n";

    C_a=rho_mat_alpha.t()*rho_mat_alpha;
    // cout << C_a;
    // cout << "\n\n";
    C_b=rho_mat_beta.t()*rho_mat_beta;
    // cout << C_b;
    // cout << "\n\n";
    P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
    P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

    P_t=P_a_new.col(0)+P_b_new.col(0);
    
    Pmat_a=P_a_new(0);
    P_aa_tot=sum(P_a_new.col(0));
    P_bb_tot=sum(P_b_new.col(0));
    Pmat_b=P_b_new(0);
    
    G_a={{g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_a,gamma_AB)},{g_mat_off(Pmat_a,gamma_AB),g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB)}};
    G_b={{g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_b,gamma_AB)},{g_mat_off(Pmat_b,gamma_AB),g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB)}};

    F_a=h_mu_nu+G_a; 
    F_b=h_mu_nu+G_b;

    energy_cdno2_new=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;
    cout << energy_cdno2_new;
    cout << "\n\n";

    cout << V_nuc;
    cout << "\n\n";

    //next iter
    fock_mat_a=F_a;
    fock_mat_b=F_b;

    eig_sym(ep_a,rho_mat_alpha,fock_mat_a);
    cout << ep_a;
    cout << "\n\n";
    cout<<rho_mat_alpha;
    cout << "\n\n";

    eig_sym(ep_b,rho_mat_beta,fock_mat_b);
    cout << ep_b;
    cout << "\n\n";
    cout<<rho_mat_beta;
    cout << "\n\n";

    C_a=rho_mat_alpha.t()*rho_mat_alpha;
    // cout << C_a;
    // cout << "\n\n";
    C_b=rho_mat_beta.t()*rho_mat_beta;
    // cout << C_b;
    // cout << "\n\n";
    P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
    P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

    P_t=P_a_new.col(0)+P_b_new.col(0);
    
    Pmat_a=P_a_new(0);
    P_aa_tot=sum(P_a_new.col(0));
    P_bb_tot=sum(P_b_new.col(0));
    Pmat_b=P_b_new(0);
    
    G_a={{g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_a,gamma_AB)},{g_mat_off(Pmat_a,gamma_AB),g_mat(P_aa_tot,Pmat_a,P_bb_tot,gamma_AA,gamma_AB)}};
    G_b={{g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB),g_mat_off(Pmat_b,gamma_AB)},{g_mat_off(Pmat_b,gamma_AB),g_mat(P_aa_tot,Pmat_b,P_bb_tot,gamma_AA,gamma_AB)}};

    F_a=h_mu_nu+G_a; 
    F_b=h_mu_nu+G_b;

    energy_cdno2_new=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;
    cout << energy_cdno2_new;
    cout << "\n\n";

    // cout << V_nuc;
    // cout << "\n\n";

  }

  return energy_cdno2_new;

}




int main(int argc, char* argv[])
{
  Shell sh1,sh2,sh3;
  string fname_1=argv[1];
  ReadShellparameter(sh1,sh2,sh3,fname_1);
  Shell shell_arr[3]={sh1,sh2,sh3};

  Shell sh4,sh5,sh6;
  string fname_2=argv[2];
  ReadShellparameter(sh4,sh5,sh6,fname_2);
  Shell shell_arr_2[3]={sh4,sh5,sh6};



  arma::mat gamma_mat={{gamma_func(shell_arr,shell_arr),gamma_func(shell_arr,shell_arr_2)},{gamma_func(shell_arr_2,shell_arr),gamma_func(shell_arr_2,shell_arr_2)}};
  // gamma_mat.print();

  double s_mu_nu=sum_func_opt(shell_arr,shell_arr_2);
  double gamma_AA=gamma_func(shell_arr,shell_arr);
  // cout << gamma_AA;
  // cout << "\n\n";
  double gamma_AB=gamma_func(shell_arr_2,shell_arr);
  // cout << gamma_AB;
  // cout << "\n\n";

  arma::mat h_core_mat={{h_core_diag(7.176,1,1,gamma_AA,gamma_AB),h_core_diag_off(-9,-9,s_mu_nu)},{h_core_diag_off(-9,-9,s_mu_nu),h_core_diag(7.176,1,1,gamma_AA,gamma_AB)}};
  // h_core_mat.print();

  arma::mat Pa=scf(h_core_mat);
  arma::mat Pa_tot=Pa.col(0)+Pa.col(0);
  // Pa_tot.print();

  arma::mat g_mat_a={{g_mat(1,0.5,1,gamma_AA,gamma_AB),g_mat_off(0.5,gamma_AB)},{g_mat_off(0.5,gamma_AB),g_mat(1,0.5,1,gamma_AA,gamma_AB)}};
  // g_mat_a.print();

  arma::mat result=g_mat_a+h_core_mat;
  // result.print();


  double v_nuc=scf_real_v_nuc(1,1,shell_arr,shell_arr_2);
  // cout << v_nuc;
  arma::mat g_mu_nu(2, 2,arma::fill::zeros);
  // cout << energy_cdno2_scf(v_nuc,gamma_AA,gamma_AB,h_core_mat,g_mu_nu);
  cout << scf_real(v_nuc,gamma_AA,gamma_AB,h_core_mat,g_mu_nu);



  
  

  return EXIT_SUCCESS;
}