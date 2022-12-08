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


double gamma_func_y(Shell* arr1,Shell* arr2)
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

          arma::vec zer_ab=pt1*(pt2-pt3)*r_ab;


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
            tot+=d_prime_1*d_prime_2*d_prime_3*d_prime_4*zer_ab(1);
          }



          

        }
      }
    }
  }
  return tot*27.2114;
}




void create_ov_mat_y(arma::mat &overlap, Shell** arr1)
{
    double row_number=size(overlap)[0];
    double col_number= size(overlap)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {
          overlap(i,j)=sum_func_opt_y(arr1[i],arr1[j]);

        }
    }

}


void create_gamma_mat_y(arma::mat &gamma, Shell** arr1)
{
    double row_number=size(gamma)[0];
    double col_number= size(gamma)[1];
    for (int i = 0; i < row_number; i++) 
    {
        for (int j = 0; j < col_number; j++) 
        {

          
          gamma(i,j)=gamma_func_y(arr1[i],arr1[j]);
          
        }
    }
}



//vnuc mat x
double vnuc_y(Shell* arr1,Shell* arr2)
{
  arma::vec Ra=arr1[0].get_R0();
  arma::vec Rb=arr2[0].get_R0();
  arma::vec r_ab=Ra-Rb;    
  arma::vec r_ab_=pow(r_ab,2);
  double R_ab_dist=sqrt(r_ab_(0)+r_ab_(1)+r_ab_(2));
  arma::vec V_nuc=(-(arr1[0].get_elem_num()*arr2[0].get_elem_num()))*r_ab/pow(R_ab_dist,3);

  return V_nuc(1)*27.2114;
}


void v_nuc_y(arma::mat &v_nuc,Shell**arr1)
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
        v_nuc(i,j)=vnuc_y(arr1[i],arr1[j]);


      }


    }

  }

}

