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

        }
    }

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

void create_gamma_diag(arma::mat &gamma, Shell** arr1,Shell** arr2)
{
    double row_number=size(gamma)[0];
    double col_number= size(gamma)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0]))
          {
            gamma(i,j)=gamma_func(arr1[i],arr2[j]);
          }
          else
          {
            gamma(i,j)=0;
          }

        }
    }
}





void hcore_sum_not(arma::mat &gamma,Shell** arr1,Shell** arr2)
{
  double row_number=size(gamma)[0];
  double col_number=size(gamma)[1];

  for(int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0]))
          {
            // gamma(i,j)=gamma_func(arr1[i],arr2[j])*arr2[j][0].get_elem_num();
            gamma(i,j)=0;
          }

          else
          {
            gamma(i,j)=gamma_func(arr1[i],arr2[j])*arr2[j][0].get_elem_num();
          }

        }

    }

}



void gamma_ptot_not_a(arma::mat &gamma,Shell** arr1,Shell** arr2)
{
  double row_number=size(gamma)[0];
  double col_number=size(gamma)[1];

  for(int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0]))
          {
            // gamma(i,j)=gamma_func(arr1[i],arr2[j])*arr2[j][0].get_elem_num();
            gamma(i,j)=0;
          }

          else
          {
            // cout << row_number/2;
            // cout << "\n\n";
            gamma(i,j)=gamma_func(arr1[i],arr2[j])*(row_number/2);
          }

        }

    }

}

void create_h_core_mat(arma::mat &hcore,arma::mat gamma2_off,arma::mat gamma3_diag,Shell** arr1)
{
    double row_number=size(hcore)[0];
    double col_number= size(hcore)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
            if(i==j)
            {
              hcore(i,j)=arr1[i][0].get_gamma()-(arr1[i][0].get_elem_num()-0.5)*sum(gamma3_diag.row(j))-sum(gamma2_off.row(j));
              // -(arr1[j][0].get_elem_num()*gamma_func(arr1[i],arr1[j]));
            }
            else
            {
              hcore(i,j)=0.5*(arr1[i][0].get_beta()+arr1[j][0].get_beta())*sum_func_opt(arr1[i],arr1[j]);
            }

        }
    }
}


double scf_real_v_nuc(double n_atom, Shell** arr1)
{
  double V_nuc=0;
  for (int i = 0; i < n_atom; i++)
  {
    for (int j = 0; j < i; j++) 
    {
      arma::vec Ra=arr1[i][0].get_R0();
      arma::vec Rb=arr1[j][0].get_R0();
      arma::vec r_ab=Ra-Rb;    
      arma::vec r_ab_=pow(r_ab,2);
      double R_ab_dist=sqrt(r_ab_(0)+r_ab_(1)+r_ab_(2));
      V_nuc+=(arr1[i][0].get_elem_num()*arr1[j][0].get_elem_num())/R_ab_dist;

    }

  }
  return V_nuc*27.2114;
}


void g_mat(arma::mat &gmat,arma::mat gamma2_off,arma::mat gamma3_diag,arma::mat pmat,Shell** arr1)
{
  double row_number=size(gmat)[0];
  double col_number= size(gmat)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i==j)
      {
        gmat(i,j)=((row_number/2)-pmat.diag()[i])*sum(gamma3_diag.row(i))+sum(gamma2_off.row(i));
      }

      else
      {
        // gmat(i,j)=0;
        gmat(i,j)=-pmat(i,j)*sum(gamma2_off.row(i));

      }

    }

  }

}



double scf_real(int n_atoms,arma::mat h_mu_nu,arma::mat g_mu_nu,arma::mat gmat_a, arma::mat gmat_b, arma::mat gamma3_offdiag,arma::mat gamma5_diag,Shell **arr1,Shell**arr3)
{
  double energy_cdno2;
  double energy_cdno2_new;
  double tol=1e-4;

  while(abs(energy_cdno2-energy_cdno2_new)<tol)
  {

    arma::mat hcore=h_mu_nu;

    arma::mat fock_mat_a=h_mu_nu+g_mu_nu;

    arma::mat fock_mat_b=h_mu_nu+g_mu_nu;

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




    //gamma3
    hcore_sum_not(gamma3_offdiag,arr3,arr1);
    //gamma5
    gamma_ptot_not_a(gamma5_diag,arr3,arr1);
  

    g_mat(gmat_a,gamma3_offdiag,gamma5_diag,P_a_new,arr1);
    g_mat(gmat_b,gamma3_offdiag,gamma5_diag,P_b_new,arr1);

      
    arma::mat G_a=gmat_a;

    arma::mat G_b=gmat_b;


    arma::mat F_a=h_mu_nu+G_a; 
    arma::mat F_b=h_mu_nu+G_b;

    
    double V_nuc = scf_real_v_nuc(n_atoms,arr1);
    cout << V_nuc;
    cout << "\n\n";

    energy_cdno2=(0.5*ep_a(0))+(0.5*ep_b(0))+V_nuc;

    //next iter
    fock_mat_a=F_a;
    fock_mat_b=F_b;


    eig_sym(ep_a,rho_mat_alpha,fock_mat_a);
    eig_sym(ep_b,rho_mat_beta,fock_mat_b);

    C_a=rho_mat_alpha.t()*rho_mat_alpha;
    C_b=rho_mat_beta.t()*rho_mat_beta;
    P_a_new=rho_mat_alpha.col(0)*rho_mat_alpha.col(0).t();
    P_b_new=rho_mat_beta.col(0)*rho_mat_beta.col(0).t();

    P_t=P_a_new.col(0)+P_b_new.col(0);

    g_mat(G_a,gamma3_offdiag,gamma5_diag,P_a_new,arr1);
    g_mat(G_b,gamma3_offdiag,gamma5_diag,P_b_new,arr1);
    
    G_a=G_a;
    G_b=G_b;

    F_a=h_mu_nu+G_a; 
    F_b=h_mu_nu+G_b;

    arma::mat Ptotal= P_a_new+P_b_new;
    double ptot=arma::dot(Ptotal, hcore) / 2;

    cout << ep_a; 
    cout << "\n\n";

    double elec_eng=ep_a(0)+ptot;
    cout << elec_eng;
    cout << "\n\n";
    
    


    energy_cdno2_new=(0.5*(ep_a(0)+ep_a(0)))+(0.5*(ep_b(0)+ep_b(0)))+V_nuc;

  }
  return energy_cdno2_new;
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
  
//C1-2px
  Shell sh10,sh11,sh12;
  string fname_4=argv[4];
  ReadShellparameter(sh10,sh11,sh12,fname_4);
  Shell shell_arr_4[3]={sh10,sh11,sh12};

//C1-2py
  Shell sh13,sh14,sh15;
  string fname_5=argv[5];
  ReadShellparameter(sh13,sh14,sh15,fname_5);
  Shell shell_arr_5[3]={sh13,sh14,sh15};

//C1-2pz
  Shell sh16,sh17,sh18;
  string fname_6=argv[6];
  ReadShellparameter(sh16,sh17,sh18,fname_6);
  Shell shell_arr_6[3]={sh16,sh17,sh18};

//C2-2s
  Shell sh19,sh20,sh21;
  string fname_7=argv[7];
  ReadShellparameter(sh19,sh20,sh21,fname_7);
  Shell shell_arr_7[3]={sh19,sh20,sh21};


//C2-2px
  Shell sh22,sh23,sh24;
  string fname_8=argv[8];
  ReadShellparameter(sh22,sh23,sh24,fname_8);
  Shell shell_arr_8[3]={sh22,sh23,sh24};

//C2-2py
  Shell sh25,sh26,sh27;
  string fname_9=argv[9];
  ReadShellparameter(sh25,sh26,sh27,fname_9);
  Shell shell_arr_9[3]={sh25,sh26,sh27};

//C2-2pz
  Shell sh28,sh29,sh30;
  string fname_10=argv[10];
  ReadShellparameter(sh28,sh29,sh30,fname_10);
  Shell shell_arr_10[3]={sh28,sh29,sh30};


//ov
  arma::mat overlap(10,10,arma::fill::zeros);
  Shell* arr2[10]={shell_arr,shell_arr_3,shell_arr_4,shell_arr_5,shell_arr_6,shell_arr_7,shell_arr_8,shell_arr_9,shell_arr_10,shell_arr_2};
  create_ov_mat(overlap,arr2);
  overlap.print();
  cout << "\n\n";


//gamma
  arma::mat gamma(4, 4,arma::fill::zeros);
  Shell* arr1[4]={shell_arr,shell_arr_3,shell_arr_7,shell_arr_2};
  create_gamma_mat(gamma,arr1);
  gamma.print();
  cout << "\n\n";


  Shell* arr3[10]={shell_arr,shell_arr_3,shell_arr_3,shell_arr_3,shell_arr_3,shell_arr_7,shell_arr_7,shell_arr_7,shell_arr_7,shell_arr_2};


//gamma3
  arma::mat gamma3_offdiag(10,4,arma::fill::zeros);

//gamma5
  arma::mat  gamma5_diag(10,4,arma::fill::zeros);


//gamma3
  arma::mat gamma3(10,4,arma::fill::zeros);
  hcore_sum_not(gamma3,arr3,arr1);

//gamma4
  arma::mat gamma4(10,4,arma::fill::zeros);
  create_gamma_diag(gamma4,arr3,arr1);

//hmat
  arma::mat hmat(10,10,arma::fill::zeros);
  create_h_core_mat(hmat,gamma3,gamma4,arr2);


//gmat test
  arma::mat gmat(10,10,arma::fill::zeros);
  arma::mat gmat_a(10,10,arma::fill::zeros);
  arma::mat gmat_b(10,10,arma::fill::zeros);


  cout << scf_real(4,hmat,gmat,gmat_a,gmat_b, gamma3_offdiag,gamma5_diag,arr1,arr3);
  cout << "\n\n";





  return EXIT_SUCCESS;
}