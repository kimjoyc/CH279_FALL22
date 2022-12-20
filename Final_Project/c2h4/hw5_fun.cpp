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
#include "new3.cpp"


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


double sum_func_opt_x(Shell* arr1,Shell* arr2)
{
    double ov_tot=0;
    for(int i=0;i<=2;i++)
    {
        for(int k=0;k<=2;k++)
        {
            double orb_1=normalize_func(arr1[i]);
            double orb_2=normalize_func(arr2[k]);
            ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1*orb_2*deriv_overlap_x(arr1[i],arr2[k])*integral_ov_y(arr1[i],arr2[k])*integral_ov_z(arr1[i],arr2[k]);

            
        }
    }
    
    return ov_tot;    

}

double sum_func_opt_y(Shell* arr1,Shell* arr2)
{
    double ov_tot=0;
    for(int i=0;i<=2;i++)
    {
        for(int k=0;k<=2;k++)
        {
            double orb_1=normalize_func(arr1[i]);
            double orb_2=normalize_func(arr2[k]);
            ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1*orb_2*integral_ov_x(arr1[i],arr2[k])*deriv_overlap_y(arr1[i],arr2[k])*integral_ov_z(arr1[i],arr2[k]);

            
        }
    }
    
    return ov_tot;    

}

double sum_func_opt_z(Shell* arr1,Shell* arr2)
{
    double ov_tot=0;
    for(int i=0;i<=2;i++)
    {
        for(int k=0;k<=2;k++)
        {
            double orb_1=normalize_func(arr1[i]);
            double orb_2=normalize_func(arr2[k]);
            ov_tot+=arr1[i].get_d()*arr2[k].get_d()*orb_1*orb_2*integral_ov_x(arr1[i],arr2[k])*integral_ov_y(arr1[i],arr2[k])*deriv_overlap_z(arr1[i],arr2[k]);  
        }
    }
    
    return ov_tot;    

}



double gamma_func_x(Shell* arr1,Shell* arr2)
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
          
          double T=V_squared*pow(R_ab_dist,2);

          double U=Ua*Ub;
          double pt1=(U/pow(R_ab_dist,2));
          double pt2=2*sqrt(V_squared)*exp(-T)/sqrt(M_PI);
          double pt3=erf(sqrt(T))/R_ab_dist;
          // arma::vec zer_AB=(U/pow(R_ab_dist,2))*((2*sqrt(V_squared)/sqrt(M_PI))*exp(T)-erf(sqrt(T))/R_ab_dist)*r_ab;

          arma::vec zer_ab=pt1*(pt2-pt3)*r_ab;
          // double zer_AB_x=zer_AB(0);
          // double zer_AB_y=zer_AB(1);
          // double zer_AB_z=zer_AB(2);


          double d_prime_1=arr1[i].get_d()*orb_1;
          double d_prime_2=arr1[k].get_d()*orb_2;
          double d_prime_3=arr2[l].get_d()*orb_3;
          double d_prime_4=arr2[m].get_d()*orb_4;


          if(R_ab_dist==0)
          {
            tot+=0;
          }
          else
          {
            tot+=d_prime_1*d_prime_2*d_prime_3*d_prime_4*zer_ab(0);
          }



          

        }
      }
    }
  }
  return tot*27.2114;
}



//ov x 
void create_ov_mat_x(arma::mat &overlap, Shell** arr1)
{
    double row_number=size(overlap)[0];
    double col_number= size(overlap)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          overlap(i,j)=sum_func_opt_x(arr1[i],arr1[j]);

        }
    }

}


//gamma x
void create_gamma_mat_x(arma::mat &gamma, Shell** arr1)
{
    double row_number=size(gamma)[0];
    double col_number= size(gamma)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {

          
          gamma(i,j)=gamma_func_x(arr1[i],arr1[j]);
          
        }
    }
}



//off diag of hcore for x
void hcore_sum_not_x(arma::mat &gamma,Shell** arr1,Shell** arr2)
{
  double row_number=size(gamma)[0];
  double col_number=size(gamma)[1];

  for(int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0]))
          {
            gamma(i,j)=0;
          }

          else
          {
            gamma(i,j)=gamma_func_x(arr1[i],arr2[j])*arr2[j][0].get_elem_num();
          }

        }

    }
}


//hcore x
void create_h_core_mat_x(arma::mat &hcore,arma::mat gamma2_off,Shell** arr1)
{
    double row_number=size(hcore)[0];
    double col_number= size(hcore)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
            if(i==j)
            {
              hcore(i,j)=-sum(gamma2_off.row(j));
            }
            else
            {
              hcore(i,j)=0.5*(arr1[i][0].get_beta()+arr1[j][0].get_beta())*sum_func_opt_x(arr1[i],arr1[j]);
            }

        }
    }
}

//density mat total not atom A for x func
void gamma_ptot_not_a_x(arma::mat&gamma,arma::mat pt,Shell** arr1,Shell**arr2)
{
  double row_number=size(gamma)[0];
  double col_number=size(gamma)[1];

  for(int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          if((arr1[i][0].get_elem_num()==arr2[j][0].get_elem_num())&&(arr1[i][0].get_R0()[0]==arr2[j][0].get_R0()[0]))
          {
            gamma(i,j)=0;
          }

          else
          {
            gamma(i,j)=(sum(pt.row(j))-arr1[j][0].get_elem_num())*gamma_func_x(arr1[i],arr2[j]);

          }

        }

    }

}


//fock mat deriv in x 
void fock_mat_x(arma::mat &gmat,arma::mat pt,arma::mat gamma2_off,arma::mat gamma3_diag,arma::mat pmat,Shell** arr1,Shell** arr2)
{
  double row_number=size(gmat)[0];
  double col_number= size(gmat)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i==j)
      {
        gmat(i,j)=sum(gamma2_off.row(i));

      }
      else
      {
        gmat(i,j)=gamma_func_x(arr2[i],arr2[j]);

      }

    }

  }

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

      else{
        x_mu_nu(i,j)=0;


      }
      // x_mu_nu(i,j)=pmat(i,j)*(arr1[i][0].get_beta()+arr1[j][0].get_beta());
      // x_mu_nu(i,j)=(arr1[i][0].get_beta()+arr1[j][0].get_beta());

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

      // y_a_b_off_diag(i,j)=x_mu_nu(i,j)*overlap(i,j);



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

      // y_a_b_off_diag(i,j)=x_mu_nu(i,j)*overlap(i,j);

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




//vnuc mat x
double vnuc_x(Shell* arr1,Shell* arr2)
{
  arma::vec Ra=arr1[0].get_R0();
  arma::vec Rb=arr2[0].get_R0();
  arma::vec r_ab=Ra-Rb;    
  arma::vec r_ab_=pow(r_ab,2);
  double R_ab_dist=sqrt(r_ab_(0)+r_ab_(1)+r_ab_(2));
  arma::vec V_nuc=(-(arr1[0].get_elem_num()*arr2[0].get_elem_num()))*r_ab/pow(R_ab_dist,3);

  return V_nuc(0)*27.2114;
}


void v_nuc_x(arma::mat &v_nuc,Shell**arr1)
{
  double row_number=size(v_nuc)[0];
  double col_number=size(v_nuc)[1];
  for (int i = 0; i < row_number; i++) 
  {
    for (int j = 0; j < col_number; j++) 
    {
      if(i==j)
      {
        v_nuc(i,j)=0;
      }

      else
      {
        v_nuc(i,j)=vnuc_x(arr1[i],arr1[j]);


      }


    }

  }

}


