
#include <iostream>
#include <armadillo>
#include <stdexcept>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <sstream>
#include <cassert>
#include "../util.h"
#include <stdlib.h>
#include <stdexcept>
#include <armadillo>
#include <stdio.h>
using namespace std;


#include "../hw5_fun_z.cpp"


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


//ov - x 
  arma::mat overlap_x(2,2,arma::fill::zeros);
  Shell* arr1[4]={shell_arr,shell_arr_2,shell_arr_2,shell_arr};
  create_ov_mat_z(overlap_x,arr1);
  overlap_x.print();
  cout << "\n\n";


//gamma 
  arma::mat gamma(2, 2,arma::fill::zeros);
  create_gamma_mat(gamma,arr1);
  gamma.print();
  cout << "\n\n";

//gamma -x
  arma::mat gamma_x(2, 2,arma::fill::zeros);
  create_gamma_mat_z(gamma_x,arr1);
  gamma_x.print();
  cout << "\n\n";


//gamma3
  arma::mat gamma3(2,2,arma::fill::zeros);
  hcore_sum_not(gamma3,arr1,arr1);
gamma3.print();
  cout << "\n\n";

//gamma4
  arma::mat gamma4(2,2,arma::fill::zeros);
  create_gamma_diag(gamma4,arr1,arr1);
  gamma4.print();
  cout << "\n\n";


//hmat
  arma::mat hmat(2,2,arma::fill::zeros);
  create_h_core_mat(hmat,gamma3,gamma4,arr1);
  hmat.print();
  cout << "\n\n";


  arma::mat g_mu_nu_(2,2,arma::fill::zeros);
  arma::mat P_a_new = scf_real_density(2,1,1,-40.569150,hmat,g_mu_nu_,arr1,arr1,arr1);
  cout << P_a_new;
  cout << "\n\n";

  arma::mat P_b_new = scf_real_density(2,1,1,-40.569150,hmat,g_mu_nu_,arr1,arr1,arr1);
  cout << P_b_new;
  cout << "\n\n";


  arma::mat g_mu_nu(2,2,arma::fill::zeros);
  arma::mat p_tot = P_a_new+P_b_new;


  arma::mat pt(2,2,arma::fill::zeros);
  pt=create_pt(pt,p_tot,arr1,arr1);
  cout << pt;
  cout << "\n\n";
  cout << sum(pt.row(0));
  cout << "\n\n";
  cout << sum(pt.row(1));
  cout << "\n\n";


  arma::mat pt_(2,2,arma::fill::zeros);
  pt_=create_pt2(pt_,p_tot,arr1);
  cout << pt_;
  cout << "\n\n";

  arma::mat vnucx(2,2,arma::fill::zeros);
  v_nuc_z(vnucx,arr1);
  vnucx.print();
  cout << "\n\n";

  cout << sum(vnucx.row(0));
  cout << "\n\n";
  cout << sum(vnucx.row(1));
  cout << "\n\n";

//yab_off
    arma::mat y_a_b_off_diag_(2,2,arma::fill::zeros);
    y_a_b_off_diag(y_a_b_off_diag_,P_a_new,P_b_new,arr1);
    y_a_b_off_diag_.print();
    cout << "\n\n";



// //make func ??? sum A, B
arma::mat sum_yab= {{0,y_a_b_off_diag_.row(0)[1]},{y_a_b_off_diag_.row(1)[0],0}};
sum_yab.print();
  cout << "\n\n";

//yab wout sum
  arma::mat yab_2(2,2,arma::fill::zeros);
  y_a_b_(yab_2,pt,P_a_new,sum_yab,arr1);
  yab_2.print();
  cout << "\n\n";


//yab
  arma::mat yab(2,2,arma::fill::zeros);
  y_a_b(yab,pt,P_a_new,sum_yab,arr1);
  yab.print();
  cout << "\n\n";

// yab*gamma_x
  arma::mat yab_gamma(2,2,arma::fill::zeros);
  y_a_bx_gamma(yab_gamma,gamma_x,yab);
  yab_gamma.print();
  cout << "\n\n";
  cout << sum(yab_gamma.row(0));
  cout << "\n\n";
  cout << sum(yab_gamma.row(1));
  cout << "\n\n";

    arma::mat xmu(2,2,arma::fill::zeros);
    x_mu_nu(xmu,p_tot,arr1);
    cout << "\n\n";

    xmu.print();
    cout << "\n\n";


    arma::mat xmu_fix_(2,2,arma::fill::zeros);
    x_mu_ov_(xmu_fix_,xmu,overlap_x,arr1);

    xmu_fix_.print();
    cout << "\n\n";

 
    arma::mat xmu_fix(2,2,arma::fill::zeros);
    x_mu_ov(xmu_fix,xmu,overlap_x,arr1);

    xmu_fix.print();
    cout << "\n\n";

    arma::mat xsum4(2,2,arma::fill::zeros);
    xmu_sum_4(xsum4,xmu_fix,arr1,arr1,arr1);
    xsum4.print();
    cout << "\n\n";

    cout << sum(xsum4.row(0));
    cout << "\n\n";
    cout << sum(xsum4.row(1));
    cout << "\n\n";


    arma::mat vnuc_grad={sum(vnucx.row(0)), sum(vnucx.row(1))};
    vnuc_grad.print();
    cout << "\n\n";

    arma::mat elec_grad={sum(xsum4.row(0))+sum(yab_gamma.row(0)),sum(yab_gamma.row(1))+sum(xsum4.row(1))};
    elec_grad.print();
    cout << "\n\n";



    cout << elec_grad+vnuc_grad;
        cout << "\n\n";







}
