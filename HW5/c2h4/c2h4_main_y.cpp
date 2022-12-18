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
#include "../hw5_fun_y.cpp"
int main(int argc, char* argv[])
{
//H1 
  Shell sh1,sh2,sh3;
  string fname_1=argv[1];
  ReadShellparameter(sh1,sh2,sh3,fname_1);
  Shell shell_arr[3]={sh1,sh2,sh3};
//H2
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

//H3
  Shell sh31,sh32,sh33;
  string fname_11=argv[11];
  ReadShellparameter(sh31,sh32,sh33,fname_11);
  Shell shell_arr_11[3]={sh31,sh32,sh33};


//H4
  Shell sh34,sh35,sh36;
  string fname_12=argv[12];
  ReadShellparameter(sh34,sh35,sh36,fname_12);
  Shell shell_arr_12[3]={sh34,sh35,sh36};


//ov - x 
  arma::mat overlap_x(12,12,arma::fill::zeros);
  Shell* arr2[12]={shell_arr,shell_arr_3,shell_arr_4,shell_arr_5,shell_arr_6,shell_arr_7,shell_arr_8,shell_arr_9,shell_arr_10,shell_arr_2,shell_arr_11,shell_arr_12};
  create_ov_mat_y(overlap_x,arr2);
  overlap_x.print();
  cout << "\n\n";


//gamma 
  arma::mat gamma(6, 6,arma::fill::zeros);
  Shell* arr1[6]={shell_arr,shell_arr_3,shell_arr_7,shell_arr_2,shell_arr_11,shell_arr_12};
  create_gamma_mat(gamma,arr1);
  gamma.print();
  cout << "\n\n";

//gamma -x
  arma::mat gamma_x(6, 6,arma::fill::zeros);
  create_gamma_mat_y(gamma_x,arr1);
  gamma_x.print();
  cout << "\n\n";

  Shell* arr3[12]={shell_arr,shell_arr_3,shell_arr_3,shell_arr_3,shell_arr_3,shell_arr_7,shell_arr_7,shell_arr_7,shell_arr_7,shell_arr_2,shell_arr_11,shell_arr_12};

//gamma3
  arma::mat gamma3(12,6,arma::fill::zeros);
  hcore_sum_not(gamma3,arr3,arr1);

//gamma4
  arma::mat gamma4(12,6,arma::fill::zeros);
  create_gamma_diag(gamma4,arr3,arr1);
//hmat
  arma::mat hmat(12,12,arma::fill::zeros);
  create_h_core_mat(hmat,gamma3,gamma4,arr2);
  hmat.print();
  cout << "\n\n";


  arma::mat g_mu_nu_(12,12,arma::fill::zeros);
  arma::mat P_a_new = scf_real_density(6,6,6,-477.150973,hmat,g_mu_nu_,arr1,arr3,arr2);
  cout << P_a_new;
  cout << "\n\n";

  arma::mat P_b_new = scf_real_density(6,6,6,-477.150973,hmat,g_mu_nu_,arr1,arr3,arr2);
  cout << P_b_new;
  cout << "\n\n";

  arma::mat p_tot = P_a_new+P_b_new;

  cout << p_tot;
  cout << "\n\n";


  arma::mat pt(6,12,arma::fill::zeros);
  pt=create_pt(pt,p_tot,arr1,arr3);
  cout << pt;
  cout << "\n\n";

  cout << sum(pt.row(0));
  cout << "\n\n";
  cout << sum(pt.row(1));
  cout << "\n\n";
  cout << sum(pt.row(2));
  cout << "\n\n";
  cout << sum(pt.row(3));
  cout << "\n\n";
  cout << sum(pt.row(4));
  cout << "\n\n";
  cout << sum(pt.row(5));
  cout << "\n\n";

  arma::mat pt_(12,12,arma::fill::zeros);
  pt_=create_pt2(pt_,p_tot,arr3);
  cout << pt_;
  cout << "\n\n";

  arma::mat vnucx(6,6,arma::fill::zeros);
  v_nuc_y(vnucx,arr1);
  vnucx.print();
  cout << "\n\n";

  cout << sum(vnucx.row(0));
  cout << "\n\n";
  cout << sum(vnucx.row(1));
  cout << "\n\n";
  cout << sum(vnucx.row(2));
  cout << "\n\n";
  cout << sum(vnucx.row(3));
  cout << "\n\n";
  cout << sum(vnucx.row(4));
  cout << "\n\n";
  cout << sum(vnucx.row(5));
  cout << "\n\n";


//yab_off
  arma::mat y_a_b_off_diag_(12,12,arma::fill::zeros);
  y_a_b_off_diag(y_a_b_off_diag_,P_a_new,P_b_new,arr2);
  y_a_b_off_diag_.print();
  cout << "\n\n";



//make func ??? sum A, B
double c1h1=y_a_b_off_diag_.row(1)[0]+y_a_b_off_diag_.row(2)[0]+y_a_b_off_diag_.row(3)[0]+y_a_b_off_diag_.row(4)[0];
double c1c2_2s=y_a_b_off_diag_.row(4)[5]+y_a_b_off_diag_.row(4)[6]+y_a_b_off_diag_.row(4)[7]+y_a_b_off_diag_.row(4)[8];
double c1c2_2p1=y_a_b_off_diag_.row(1)[5]+y_a_b_off_diag_.row(1)[6]+y_a_b_off_diag_.row(1)[7]+y_a_b_off_diag_.row(1)[8];
double c1c2_2p2=y_a_b_off_diag_.row(2)[5]+y_a_b_off_diag_.row(2)[6]+y_a_b_off_diag_.row(2)[7]+y_a_b_off_diag_.row(2)[8];
double c1c2_2p3=y_a_b_off_diag_.row(3)[5]+y_a_b_off_diag_.row(3)[6]+y_a_b_off_diag_.row(3)[7]+y_a_b_off_diag_.row(3)[8];
double c1c2_tot=c1c2_2s+c1c2_2p1+c1c2_2p2+c1c2_2p3;
double c1h2=y_a_b_off_diag_.row(1)[9]+y_a_b_off_diag_.row(2)[9]+y_a_b_off_diag_.row(3)[9]+y_a_b_off_diag_.row(4)[9];

double c1h3=y_a_b_off_diag_.row(1)[10]+y_a_b_off_diag_.row(2)[10]+y_a_b_off_diag_.row(3)[10]+y_a_b_off_diag_.row(4)[10];
double c1h4=y_a_b_off_diag_.row(1)[11]+y_a_b_off_diag_.row(2)[11]+y_a_b_off_diag_.row(3)[11]+y_a_b_off_diag_.row(4)[11];



double c2h1=y_a_b_off_diag_.row(5)[0]+y_a_b_off_diag_.row(6)[0]+y_a_b_off_diag_.row(7)[0]+y_a_b_off_diag_.row(8)[0];
double c2c1_2s=y_a_b_off_diag_.row(5)[1]+y_a_b_off_diag_.row(5)[2]+y_a_b_off_diag_.row(5)[3]+y_a_b_off_diag_.row(5)[4];
double c2c1_2p1=y_a_b_off_diag_.row(6)[1]+y_a_b_off_diag_.row(6)[2]+y_a_b_off_diag_.row(6)[3]+y_a_b_off_diag_.row(6)[4];
double c2c1_2p2=y_a_b_off_diag_.row(7)[1]+y_a_b_off_diag_.row(7)[2]+y_a_b_off_diag_.row(7)[3]+y_a_b_off_diag_.row(7)[4];
double c2c1_2p3=y_a_b_off_diag_.row(8)[1]+y_a_b_off_diag_.row(8)[2]+y_a_b_off_diag_.row(8)[3]+y_a_b_off_diag_.row(8)[4];
double c2c1_tot=c2c1_2s+c2c1_2p1+c2c1_2p2+c2c1_2p3;
double c2h2=y_a_b_off_diag_.row(5)[9]+y_a_b_off_diag_.row(6)[9]+y_a_b_off_diag_.row(7)[9]+y_a_b_off_diag_.row(8)[9];

double c2h3=y_a_b_off_diag_.row(5)[10]+y_a_b_off_diag_.row(6)[10]+y_a_b_off_diag_.row(7)[10]+y_a_b_off_diag_.row(8)[10];
double c2h4=y_a_b_off_diag_.row(5)[11]+y_a_b_off_diag_.row(6)[11]+y_a_b_off_diag_.row(7)[11]+y_a_b_off_diag_.row(8)[11];




double h1c1=y_a_b_off_diag_.row(0)[1]+y_a_b_off_diag_.row(0)[2]+y_a_b_off_diag_.row(0)[3]+y_a_b_off_diag_.row(0)[4];
double h1c2=y_a_b_off_diag_.row(0)[5]+y_a_b_off_diag_.row(0)[6]+y_a_b_off_diag_.row(0)[7]+y_a_b_off_diag_.row(0)[8];

double h1h2=y_a_b_off_diag_.row(0)[9];
double h1h3=y_a_b_off_diag_.row(0)[10];
double h1h4=y_a_b_off_diag_.row(0)[11];



double h2h1=y_a_b_off_diag_.row(9)[0];
double h2h3=y_a_b_off_diag_.row(9)[10];
double h2h4=y_a_b_off_diag_.row(9)[11];


double h2c1=y_a_b_off_diag_.row(9)[1]+y_a_b_off_diag_.row(9)[2]+y_a_b_off_diag_.row(9)[3]+y_a_b_off_diag_.row(9)[4];
double h2c2=y_a_b_off_diag_.row(9)[5]+y_a_b_off_diag_.row(9)[6]+y_a_b_off_diag_.row(9)[7]+y_a_b_off_diag_.row(9)[8];

double h3h1=y_a_b_off_diag_.row(10)[0];
double h3h2=y_a_b_off_diag_.row(10)[9];
double h3h4=y_a_b_off_diag_.row(10)[11];


double h3c1=y_a_b_off_diag_.row(10)[1]+y_a_b_off_diag_.row(10)[2]+y_a_b_off_diag_.row(10)[3]+y_a_b_off_diag_.row(10)[4];
double h3c2=y_a_b_off_diag_.row(10)[5]+y_a_b_off_diag_.row(10)[6]+y_a_b_off_diag_.row(10)[7]+y_a_b_off_diag_.row(10)[8];


double h4h1=y_a_b_off_diag_.row(11)[0];
double h4h2=y_a_b_off_diag_.row(11)[9];
double h4h3=y_a_b_off_diag_.row(11)[10];


double h4c1=y_a_b_off_diag_.row(11)[1]+y_a_b_off_diag_.row(11)[2]+y_a_b_off_diag_.row(11)[3]+y_a_b_off_diag_.row(11)[4];
double h4c2=y_a_b_off_diag_.row(11)[5]+y_a_b_off_diag_.row(11)[6]+y_a_b_off_diag_.row(11)[7]+y_a_b_off_diag_.row(11)[8];



arma::mat sum_yab=
{
  {0
  ,y_a_b_off_diag_.row(0)[1]+y_a_b_off_diag_.row(0)[2]+y_a_b_off_diag_.row(0)[3]+y_a_b_off_diag_.row(0)[4]
  ,y_a_b_off_diag_.row(0)[5]+y_a_b_off_diag_.row(0)[6]+y_a_b_off_diag_.row(0)[7]+y_a_b_off_diag_.row(0)[8]
  ,y_a_b_off_diag_.row(0)[9]
  ,h1h3
  ,h1h4},
  {c1h1,0,c1c2_tot,c1h2,c1h3,c1h4},
  {c2h1,c2c1_tot,0,c2h2,c2h3,c2h4},
  {y_a_b_off_diag_.row(9)[0],
  y_a_b_off_diag_.row(9)[1]+y_a_b_off_diag_.row(9)[2]+y_a_b_off_diag_.row(9)[3]+y_a_b_off_diag_.row(9)[4],
  y_a_b_off_diag_.row(9)[5]+y_a_b_off_diag_.row(9)[6]+y_a_b_off_diag_.row(9)[7]+y_a_b_off_diag_.row(9)[8],
  0
  ,h2h3
  ,h2h4},

  {h3h1,h3c1,h3c2,h3h2,0,h3h4},

  {h4h1,h4c1,h4c2,h4h2,h4h3,0}

};


sum_yab.print();
cout << "\n\n";

//yab wout sum
  arma::mat yab_2(6,6,arma::fill::zeros);
  y_a_b_(yab_2,pt,P_a_new,sum_yab,arr1);
  yab_2.print();
  cout << "\n\n";


//yab
  arma::mat yab(6,6,arma::fill::zeros);
  y_a_b(yab,pt,P_a_new,sum_yab,arr1);
  yab.print();
  cout << "\n\n";

// yab*gamma_x
  arma::mat yab_gamma(6,6,arma::fill::zeros);
  y_a_bx_gamma(yab_gamma,gamma_x,yab);
  yab_gamma.print();
  cout << "\n\n";

  cout << sum(yab_gamma.row(0));
  cout << "\n\n";
  cout << sum(yab_gamma.row(1));
  cout << "\n\n";
  cout << sum(yab_gamma.row(2));
  cout << "\n\n";
  cout << sum(yab_gamma.row(3));
  cout << "\n\n";
  cout << sum(yab_gamma.row(4));
  cout << "\n\n";
  cout << sum(yab_gamma.row(5));
  cout << "\n\n";


// void x_mu_nu(arma::mat &x_mu_nu,arma::mat pmat,Shell** arr1)
    arma::mat xmu(12,12,arma::fill::zeros);
    x_mu_nu(xmu,p_tot,arr2);
    xmu.print();
    cout << "\n\n";


   // void x_mu(arma::mat &y_a_b_off_diag,arma::mat pmat,arma::mat x_mu_nu,arma::mat overlap,Shell** arr1)
    arma::mat xmu_fix_(12,12,arma::fill::zeros);
    x_mu_ov_(xmu_fix_,xmu,overlap_x,arr2);
    xmu_fix_.print();
    cout << "\n\n";

 
// void x_mu(arma::mat &y_a_b_off_diag,arma::mat pmat,arma::mat x_mu_nu,arma::mat overlap,Shell** arr1)
    arma::mat xmu_fix(12,12,arma::fill::zeros);
    x_mu_ov(xmu_fix,xmu,overlap_x,arr2);
    xmu_fix.print();
    cout << "\n\n";



// void xmu_sum_4(arma::mat&pt, arma::mat pt_tot,Shell**arr1,Shell**arr2,Shell**arr3)
    arma::mat xsum4(6,12,arma::fill::zeros);
    xmu_sum_4(xsum4,xmu_fix,arr1,arr3,arr2);
    xsum4.print();
    cout << "\n\n";

    cout << sum(xsum4.row(0));
    cout << "\n\n";
    cout << sum(xsum4.row(1));
    cout << "\n\n";
    cout << sum(xsum4.row(2));
    cout << "\n\n";
    cout << sum(xsum4.row(3));
    cout << "\n\n";
    cout << sum(xsum4.row(4));
    cout << "\n\n";
    cout << sum(xsum4.row(5));
    cout << "\n\n";



  arma::mat vnuc_grad={sum(vnucx.row(0)), sum(vnucx.row(1)), sum(vnucx.row(2)), sum(vnucx.row(3)),sum(vnucx.row(4)),sum(vnucx.row(5))};
  vnuc_grad.print();
  cout << "\n\n"; 


  double elec_grad_1=sum(yab_gamma.row(0))+sum(xsum4.row(0));
  double elec_grad_2=sum(yab_gamma.row(1))+sum(xsum4.row(1));
  double elec_grad_3=sum(yab_gamma.row(2))+sum(xsum4.row(2));
  double elec_grad_4=sum(yab_gamma.row(3))+sum(xsum4.row(3));

  double elec_grad_5=sum(yab_gamma.row(4))+sum(xsum4.row(4));
  double elec_grad_6=sum(yab_gamma.row(5))+sum(xsum4.row(5));





  arma::mat elec_grad={sum(yab_gamma.row(0))+sum(xsum4.row(0))
  ,sum(yab_gamma.row(1))+sum(xsum4.row(1))
  ,sum(yab_gamma.row(2))+sum(xsum4.row(2))
  ,sum(yab_gamma.row(3))+sum(xsum4.row(3))
  ,elec_grad_5,elec_grad_6};

  elec_grad.print();
  cout << "\n\n";
  
    cout << elec_grad+vnuc_grad;
    cout << "\n\n";






  return EXIT_SUCCESS;
}