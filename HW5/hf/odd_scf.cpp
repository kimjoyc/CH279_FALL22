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

//gamma diag (hcore)
void create_gamma_diag(arma::mat &gamma, Shell** arr1,Shell** arr2)
{
    double row_number=size(gamma)[0];
    double col_number= size(gamma)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0])&&(arr1[i][0].get_R0()[1]==arr2[j][0].get_R0()[1])&&(arr1[i][0].get_R0()[2]==arr2[j][0].get_R0()[2]))
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

//off diag of hcore
void hcore_sum_not(arma::mat &gamma,Shell** arr1,Shell** arr2)
{
  double row_number=size(gamma)[0];
  double col_number=size(gamma)[1];

  for(int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0])&&(arr1[i][0].get_R0()[1]==arr2[j][0].get_R0()[1])&&(arr1[i][0].get_R0()[2]==arr2[j][0].get_R0()[2]))
          {
            gamma(i,j)=0;
          }

          else
          {
            gamma(i,j)=arr2[j][0].get_elem_num()*gamma_func(arr1[i],arr2[j]);
          }

        }

    }
}

arma::mat gamma_ptot_not_a2(arma::mat gamma,arma::mat pt,Shell** arr1,Shell**arr2)
{
  double row_number=size(gamma)[0];
  double col_number=size(gamma)[1];

  for(int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0])&&(arr1[i][0].get_R0()[1]==arr2[j][0].get_R0()[1])&&(arr1[i][0].get_R0()[2]==arr2[j][0].get_R0()[2]))
          {
            gamma(i,j)=0;
          }

          else
          {
            gamma(i,j)=sum(pt.row(j))*gamma_func(arr1[i],arr2[j]);

          }

        }

    }

  return gamma;

}

//pt so that can actually iterate 
arma::mat create_pt(arma::mat pt, arma::mat pt_tot,Shell**arr1,Shell**arr2)
{
    double row_number=size(pt)[0];
    double col_number= size(pt)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0])&&(arr1[i][0].get_R0()[1]==arr2[j][0].get_R0()[1])&&(arr1[i][0].get_R0()[2]==arr2[j][0].get_R0()[2]))
          {
            pt(i,j)=pt_tot.diag()[j];
          }
          else
          {
            pt(i,j)=0;
          }

        }

    }

  return pt;


}

//pt so that can actually iterate 
arma::mat create_pt2(arma::mat pt2, arma::mat pt_tot,Shell**arr1)
{
    double row_number=size(pt2)[0];
    double col_number= size(pt2)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr1[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr1[j][0].get_R0()[0])&&(arr1[i][0].get_R0()[1]==arr1[j][0].get_R0()[1])&&(arr1[i][0].get_R0()[2]==arr1[j][0].get_R0()[2]))
          {
            pt2(i,j)=pt_tot.diag()[j];
          }
          else
          {
            pt2(i,j)=0;
          }

        }

    }
    return pt2;
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


arma::mat g_mat2(arma::mat gmat,arma::mat pt,arma::mat gamma2_off,arma::mat gamma3_diag,arma::mat pmat,Shell** arr1,Shell** arr2)
{
  double row_number=size(gmat)[0];
  double col_number= size(gmat)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i==j)
      {

        gmat(i,j)=(sum(pt.row(i))-pmat(i,j))*gamma_func(arr2[i],arr2[j])+sum(gamma2_off.row(i));


      }

      else
      {
        gmat(i,j)=-pmat(i,j)*gamma_func(arr2[i],arr2[j]);

      }


    }

  }
  return gmat;


}

arma::mat scf_real_density(int n_atoms,int p,int q,double tot_eng,arma::mat h_mu_nu,arma::mat g_mu_nu,Shell**arr1,Shell**arr3,Shell**arr2)
{

  double energy_cdno2=tot_eng;
  double energy_cdno2_new;
  double energy_cdno2_init;
  double iteration=0;
  double tol=1e-4;
  arma::mat Pa_new;

    //init cond
    arma::mat hcore=h_mu_nu;
    arma::mat fock_mat_a=h_mu_nu+g_mu_nu;
    arma::mat fock_mat_b=h_mu_nu+g_mu_nu;

    arma::vec ep_a;
    arma::mat rho_mat_alpha;
    eig_sym(ep_a,rho_mat_alpha,fock_mat_a);

    arma::vec ep_b;
    arma::mat rho_mat_beta;
    eig_sym(ep_b,rho_mat_beta,fock_mat_b);

    arma::mat P_a_new=rho_mat_alpha.cols(0,p-1)*rho_mat_alpha.cols(0,p-1).t();
    arma::mat P_b_new=rho_mat_beta.cols(0,q-1)*rho_mat_beta.cols(0,q-1).t();
    arma::mat P_t=P_a_new+P_b_new;

    //gamma4
    arma::mat gamma4(p+q,n_atoms,arma::fill::zeros);
    create_gamma_diag(gamma4,arr3,arr1);

  //dynamic 
    arma::mat pt(n_atoms,p+q,arma::fill::zeros);
    pt=create_pt(pt,P_t,arr1,arr3);

    arma::mat pt_(p+q,p+q,arma::fill::zeros);
    pt_=create_pt2(pt_,P_t,arr3);

    arma::mat gamma_fun2(p+q,n_atoms,arma::fill::zeros);
    gamma_fun2=gamma_ptot_not_a2(gamma_fun2,pt,arr3,arr1);

    arma::mat gmat_check_alpha(p+q,p+q,arma::fill::zeros);
    gmat_check_alpha=g_mat2(gmat_check_alpha,pt_,gamma_fun2,gamma4,P_a_new,arr2,arr3);

    arma::mat gmat_check_beta(p+q,p+q,arma::fill::zeros);
    gmat_check_beta=g_mat2(gmat_check_beta,pt_,gamma_fun2,gamma4,P_b_new,arr2,arr3);


    arma::mat G_a=gmat_check_alpha;
    arma::mat G_b=gmat_check_beta;

    arma::mat F_a=h_mu_nu+G_a; 
    arma::mat F_b=h_mu_nu+G_b;
  

    double V_nuc = scf_real_v_nuc(n_atoms,arr1);
    double elec_eng0=0.5*trace(P_a_new*(hcore+F_a))+0.5*trace(P_b_new*(hcore+F_b));
    energy_cdno2_init=elec_eng0+V_nuc;


    //next iter
    arma::mat fock_mat_a_next=F_a;
    arma::mat fock_mat_b_next=F_b;

  while(abs(energy_cdno2-energy_cdno2_new)>=tol)
  {

    arma::vec ep_a_next;
    arma::mat rho_mat_alpha_next;
    eig_sym(ep_a_next,rho_mat_alpha_next,fock_mat_a_next);

    arma::vec ep_b_next;
    arma::mat rho_mat_beta_next;
    eig_sym(ep_b_next,rho_mat_beta_next,fock_mat_b_next);

    arma::mat P_a_new_next=rho_mat_alpha_next.cols(0,p-1)*rho_mat_alpha_next.cols(0,p-1).t();
    arma::mat P_b_new_next=rho_mat_beta_next.cols(0,q-1)*rho_mat_beta_next.cols(0,q-1).t();

    arma::mat P_t_next=P_a_new_next+P_b_new_next;

  
  //gamma4
    arma::mat gamma4_next(p+q,n_atoms,arma::fill::zeros);
    create_gamma_diag(gamma4_next,arr3,arr1);

    arma::mat pt_next(n_atoms,p+q,arma::fill::zeros);
    pt_next=create_pt(pt_next,P_t_next,arr1,arr3);

    arma::mat nextpt_(p+q,p+q,arma::fill::zeros);
    nextpt_=create_pt2(nextpt_,P_t_next,arr3);

    arma::mat gamma_fun2_next(p+q,n_atoms,arma::fill::zeros);
    gamma_fun2_next=gamma_ptot_not_a2(gamma_fun2_next,pt_next,arr3,arr1);

    arma::mat gmat_check_alpha_next(p+q,p+q,arma::fill::zeros);
    gmat_check_alpha_next=g_mat2(gmat_check_alpha_next,nextpt_,gamma_fun2_next,gamma4_next,P_a_new_next,arr2,arr3);

    arma::mat gmat_check_beta_next(p+q,p+q,arma::fill::zeros);
    gmat_check_beta_next=g_mat2(gmat_check_beta_next,nextpt_,gamma_fun2_next,gamma4_next,P_b_new_next,arr2,arr3);


    arma::mat G_a_next=gmat_check_alpha_next;
    arma::mat G_b_next=gmat_check_beta_next;

    arma::mat F_a_next=h_mu_nu+G_a_next; 
    arma::mat F_b_next=h_mu_nu+G_b_next;

    double V_nuc_next=scf_real_v_nuc(n_atoms,arr1);
    double elec_eng1=0.5*trace(P_a_new_next*(hcore+F_a_next))+0.5*trace(P_b_new_next*(hcore+F_b_next));
    energy_cdno2_new=elec_eng1+V_nuc_next;



    //next iter2
    fock_mat_a_next=F_a_next;
    fock_mat_b_next=F_b_next;

    Pa_new=P_a_new_next;

  }


  return Pa_new;


}




double integral_ov_x(Shell &sh1, Shell & sh2)
{

  arma::vec la = sh1.get_l();
  arma::vec lb = sh2.get_l();

  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();

  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();

  double overlap_int_x=Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la(0), lb(0));
  return overlap_int_x;
}


double integral_ov_y(Shell &sh1, Shell & sh2)
{

  arma::vec la = sh1.get_l();
  arma::vec lb = sh2.get_l();

  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();

  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();

  double overlap_int_y=Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la(1), lb(1));
  return overlap_int_y;
}

double integral_ov_z(Shell &sh1, Shell & sh2)
{

  arma::vec la = sh1.get_l();
  arma::vec lb = sh2.get_l();

  double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();

  arma::vec Ra = sh1.get_R0();
  arma::vec Rb = sh2.get_R0();

  double overlap_int_z=Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la(2),lb(2));
  return overlap_int_z;
}

double deriv_overlap_x(Shell &sh1, Shell &sh2)
{
    arma::vec la = sh1.get_l();
    arma::vec lb = sh2.get_l();

    double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();

    arma::vec Ra = sh1.get_R0();
    arma::vec Rb = sh2.get_R0();


    double overlap_x_der= (-la(0)*Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la(0)-1, lb(0)))+((2*alphaa)*
    Overlap_onedim(Ra(0), Rb(0), alphaa, alphab, la(0)+1, lb(0)));
    return overlap_x_der;

}


double deriv_overlap_y(Shell &sh1, Shell &sh2)
{
    arma::vec la = sh1.get_l();
    arma::vec lb = sh2.get_l();

    double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
    arma::vec Ra = sh1.get_R0();
    arma::vec Rb = sh2.get_R0();


    double overlap_y_der= -la(1)*Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la(1)-1, lb(1))+((2*alphaa)*
    Overlap_onedim(Ra(1), Rb(1), alphaa, alphab, la(1)+1, lb(1)));

    return overlap_y_der;
}

double deriv_overlap_z(Shell &sh1, Shell &sh2)
{
    arma::vec la = sh1.get_l();
    arma::vec lb = sh2.get_l();

    double alphaa = sh1.get_alpha(), alphab = sh2.get_alpha();
    arma::vec Ra = sh1.get_R0();
    arma::vec Rb = sh2.get_R0();

    double overlap_z_der= -la(2)*Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la(2)-1, lb(2))+((2*alphaa)*
    Overlap_onedim(Ra(2), Rb(2), alphaa, alphab, la(2)+1, lb(2)));

    return overlap_z_der;

}




//y_a_b gradient for gamma summation matrix offdiag
void y_a_b_off_diag(arma::mat &y_a_b_off_diag,arma::mat pmat_alpha,arma::mat pmat_beta,Shell** arr1)
{
  double row_number=size(y_a_b_off_diag)[0];
  double col_number=size(y_a_b_off_diag)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      y_a_b_off_diag(i,j)=pmat_alpha(i,j)*pmat_alpha(i,j)+pmat_beta(i,j)*pmat_beta(i,j);
    }


  }
}




//y_a_b gradient for gamma
void y_a_b_(arma::mat &y_a_b,arma::mat p_tot,arma::mat pmat,arma::mat pmat_off_diag,Shell** arr1)
{
  double row_number=size(y_a_b)[0];
  double col_number=size(y_a_b)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i==j)
      {
        y_a_b(i,j)=0;
      }

      else
      {
        y_a_b(i,j)=sum(p_tot.row(i))*sum(p_tot.row(j))-arr1[i][0].get_elem_num()*sum(p_tot.row(j))-arr1[j][0].get_elem_num()*sum(p_tot.row(i));


      }
      // -pmat_off_diag(i,j);

    }

  }

}


//y_a_b gradient for gamma
void y_a_b(arma::mat &y_a_b,arma::mat p_tot,arma::mat pmat,arma::mat pmat_off_diag,Shell** arr1)
{
  double row_number=size(y_a_b)[0];
  double col_number=size(y_a_b)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i!=j)
      {
        y_a_b(i,j)=sum(p_tot.row(i))*sum(p_tot.row(j))-arr1[i][0].get_elem_num()*sum(p_tot.row(j))-arr1[j][0].get_elem_num()*sum(p_tot.row(i))-pmat_off_diag(i,j);
      }

      else{
        y_a_b(i,j)=0;
      }

    }

  }

}



//y_a_b gradient for gamma * gamma_x
void y_a_bx_gamma(arma::mat &y_a_b_gamma,arma::mat gamma,arma::mat yab)
{
  double row_number=size(y_a_b_gamma)[0];
  double col_number=size(y_a_b_gamma)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      y_a_b_gamma(i,j)=gamma(i,j)*yab(i,j);
      
    }

  }

}


//gradient of x_mu_nu
void x_mu_nu(arma::mat &x_mu_nu,arma::mat pmat,Shell** arr1)
{
  double row_number=size(x_mu_nu)[0];
  double col_number=size(x_mu_nu)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i!=j)
      {
        x_mu_nu(i,j)=pmat(i,j)*(arr1[i][0].get_beta()+arr1[j][0].get_beta());


      }

      else
      {
        x_mu_nu(i,j)=0;

      }

    }

  }

}


//xmu gradient for gamma summation matrix offdiag
void x_mu_ov_(arma::mat &y_a_b_off_diag,arma::mat x_mu_nu,arma::mat overlap,Shell** arr1)
{
  double row_number=size(y_a_b_off_diag)[0];
  double col_number=size(y_a_b_off_diag)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if((arr1[i][0].get_elem_num()==arr1[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr1[j][0].get_R0()[0]))
      {
        y_a_b_off_diag(i,j)=0;

      }


      else
      {
        y_a_b_off_diag(i,j)=overlap(i,j);

      }




    }

  }
}


//xmu gradient for gamma summation matrix offdiag
void x_mu_ov(arma::mat &y_a_b_off_diag,arma::mat x_mu_nu,arma::mat overlap,Shell** arr1)
{
  double row_number=size(y_a_b_off_diag)[0];
  double col_number=size(y_a_b_off_diag)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if((arr1[i][0].get_elem_num()==arr1[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr1[j][0].get_R0()[0]))
      {
        y_a_b_off_diag(i,j)=0;

      }


      else
      {
        y_a_b_off_diag(i,j)=x_mu_nu(i,j)*overlap(i,j);

      }


    }

  }
}


//pt so that can actually iterate 
void xmu_sum_4(arma::mat&pt, arma::mat pt_tot,Shell**arr1,Shell**arr2,Shell**arr3)
{
    double row_number=size(pt)[0];
    double col_number= size(pt)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0]))
          {
            pt(i,j)=sum(pt_tot.row(j));
          }
          else
          {
            pt(i,j)=0;
          }


        }

    }

}





