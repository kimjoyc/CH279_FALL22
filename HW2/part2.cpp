#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>      
#include <armadillo>


using std::string;
using std::vector;
using namespace std;

class S_AB{
private:
	float L_A;
    float L_B;
	float R_A;
	float R_B;
    float alpha_;
    float beta_;

public:
	S_AB(float L_A_, float L_B_, float R_A_, float R_B_, float alpha_2,float beta_2) {
        L_A=L_A_;
        L_B=L_B_;
        R_A=R_A_;
        R_B=R_B_;
        alpha_=alpha_2;
        beta_=beta_2;
	}

    /*You need to implement the factorial and double factorial functions. This can be done as a loop, where you might most simply have separate cases for odd and even n!!. */
    float double_factorial(float n)
    {
        if (n == 0||n==1||n==-1)
            return 1;

        return n*double_factorial(n-2);

    }

    float factorial(float n) 
    {
        if ((n==0)||(n==1)){
            return 1;
        }
        else{
            return n*factorial(n-1);
        }
    }

    /*Evaluate the center of the product gaussian, RP */
    float gaussian_prod_center(){
        float R_P=(alpha_*R_A+beta_*R_B)/(alpha_+beta_);
        return R_P;
    }

    /*Evaluate the exponential prefactor of Eq. 2.9, and the associated square root.*/
    float exp_prefecator(){
        return sqrt(M_PI/(alpha_+beta_));
    }
    /* Decide how general your code should be. At a minimum, we need it to work for s shells and p shells. If you can get it to work for higher angular momentum shells as well, that would be cool, but is not essential. */
    float gaussian_omega_r(float l, float x,float X_center,float alpha){
        return pow(x-X_center,l)*exp(-alpha*pow(x-X_center,l));

    }

    float gaussian_alpha_beta_exp(){
        return exp(-(alpha_*beta_*pow(R_A-R_B,2))/(alpha_+beta_));
    }

    /*Write some code to execute the double summation for a pair of l A and lB values associ-ated with the shell pair. For instance, for a p −s shell pair, l A =0, 1 and lB =0 meaning there are 2 possible (l A , lB ) pairs, while there is only one for an s −s shell pair.*/
    float double_summation()
    {
        float total=0;
        for(int i=0; i<=L_A; i++)
        {
            for(int j=0; j<=L_B; j++)
            {
                if (i+j%2==0)
                {
                    float fact_diff_B= L_B-j;
                    float fact_tot_B = factorial(L_B)/factorial(j)*factorial(fact_diff_B);
                    float fact_diff_A= L_A-i;
                    float fact_tot_A = factorial(L_A)/factorial(i)*factorial(fact_diff_A);
                    float n = i + j - 1; 
                    float increment_tot=double_factorial(n)*pow(gaussian_prod_center()-R_A,L_A-i)*pow(gaussian_prod_center()-R_B,L_B-j)/pow(2*(alpha_+beta_),(i+j)/2);
                    float sum = gaussian_alpha_beta_exp()*exp_prefecator()*increment_tot*fact_tot_A*fact_tot_B;
                    total += sum;
                    return total;
                }
                else 
                    return total;
            }

        }
        return total;

    }

};

int main()
{
   /*Check your analytical integration code for S ABx against numerical integration from the first problem! Check that cases which should give zero do give zero. We will discuss some of these issues in the compute lab.*/
    // Analytical Integral 1.25331413731550012e+00 for s1s1

    //This shell info: R( 0.00, 0.00, 0.00), with angular momentum: 1, coefficient: 1.00
    //This shell info: R( 1.00, 1.00, 2.00), with angular momentum: 1, coefficient: 1.00 (SHELL2 )

   float L_A_s1_x=0;
   float L_B_s1_x=0;
   float R_A_s1_x=0;
   float R_B_s1_x=0;
   float alpha_s1_x=1;
   float beta_s1_x=1;
   S_AB s1s1_x(L_A_s1_x,L_B_s1_x,R_A_s1_x,R_B_s1_x,alpha_s1_x,beta_s1_x);


   //tot
   cout << "s1s1 tot " << s1s1_x.double_summation();
   cout << "\n\n";


    // Analytical Integral 0.00000000000000000e+00 for s1p1
    //x
   float L_Apx_x=1;
   float L_Bs_x=0;
   float R_Apx_x=0;
   float R_Bs_x=0;
   float alpha_px_x=1;
   float beta_s_x=1;
   S_AB s1p1x_x(L_Apx_x,L_Bs_x,R_Apx_x,R_Bs_x,alpha_px_x,beta_s_x);
   cout << "x "<< s1p1x_x.double_summation();
   cout << "\n\n";
   
   
   float L_Apx_y=0;
   float L_Bs_y=0;
   float R_Apx_y=0;
   float R_Bs_y=0;
   float alpha_px_y=1;
   float beta_s_y=1;
   S_AB s1p1x_y(L_Apx_y,L_Bs_y,R_Apx_y,R_Bs_y,alpha_px_y,beta_s_y);
   cout << "y "<< s1p1x_y.double_summation();
   cout << "\n\n";
   
   
   float L_Apx_z=0;
   float L_Bs_z=0;
   float R_Apx_z=0;
   float R_Bs_z=0;
   float alpha_px_z=1;
   float beta_s_z=1;
   S_AB s1p1x_z(L_Apx_z,L_Bs_z,R_Apx_z,R_Bs_z,alpha_px_z,beta_s_z);
   cout << "z "<< s1p1x_z.double_summation();
   cout << "\n\n";
   
   //x tot
    float s1p1_x_tot=s1p1x_x.double_summation()*s1p1x_y.double_summation()*s1p1x_z.double_summation();

   cout << "s1p1_p1:  " << s1p1x_x.double_summation()*s1p1x_y.double_summation()*s1p1x_z.double_summation();
   cout << "\n\n";

   //y
   float L_Apy_x=0;
//    float L_Bs_x=0;
   float R_Apy_x=0;
//    float R_Bs_x=0;
   float alpha_py_x=1;
//    float beta_s_x=1;
   S_AB s1p1y_x(L_Apy_x,L_Bs_x,R_Apy_x,R_Bs_x,alpha_py_x,beta_s_x);
   cout << "x "<< s1p1y_x.double_summation();
   cout << "\n\n";
   
   
   float L_Apy_y=1;
//    float L_Bs_y=0;
   float R_Apy_y=0;
//    float R_Bs_y=0;
   float alpha_py_y=1;
//    float beta_s_y=1;
   S_AB s1p1y_y(L_Apy_y,L_Bs_y,R_Apy_y,R_Bs_y,alpha_py_y,beta_s_y);
   cout << "y "<< s1p1y_y.double_summation();
   cout << "\n\n";
   
   
   float L_Apy_z=0;
//    float L_Bs_z=0;
   float R_Apy_z=0;
//    float R_Bs_z=0;
   float alpha_py_z=1;
//    float beta_s_z=1;
   S_AB s1p1y_z(L_Apy_z,L_Bs_z,R_Apy_z,R_Bs_z,alpha_py_z,beta_s_z);
   cout << "z "<< s1p1y_z.double_summation();
   cout << "\n\n";
   
   //y tot
   float s1p1_y_tot=s1p1y_x.double_summation()*s1p1y_y.double_summation()*s1p1y_z.double_summation();
   cout << "s1p1_p2:  " << s1p1y_x.double_summation()*s1p1y_y.double_summation()*s1p1y_z.double_summation();
   cout << "\n\n";



   //z
   float L_Apz_x=0;
//    float L_Bs_x=0;
   float R_Apz_x=0;
//    float R_Bs_x=0;
   float alpha_pz_x=1;
//    float beta_s_x=1;
   S_AB s1p1z_x(L_Apz_x,L_Bs_x,R_Apz_x,R_Bs_x,alpha_pz_x,beta_s_x);
   cout << "x "<< s1p1z_x.double_summation();
   cout << "\n\n";
   
   
   float L_Apz_y=0;
//    float L_Bs_y=0;
   float R_Apz_y=0;
//    float R_Bs_y=0;
   float alpha_pz_y=1;
//    float beta_s_y=1;
   S_AB s1p1z_y(L_Apz_y,L_Bs_y,R_Apz_y,R_Bs_y,alpha_pz_y,beta_s_y);
   cout << "y "<< s1p1z_y.double_summation();
   cout << "\n\n";
   
   
   float L_Apz_z=1;
//    float L_Bs_z=0;
   float R_Apz_z=0;
//    float R_Bs_z=0;
   float alpha_pz_z=1;
//    float beta_s_z=1;
   S_AB s1p1z_z(L_Apz_z,L_Bs_z,R_Apz_z,R_Bs_z,alpha_pz_z,beta_s_z);
   cout << "z "<< s1p1z_z.double_summation();
   cout << "\n\n";
   
   //z tot
   float s1p1_z_tot=s1p1z_x.double_summation()*s1p1z_y.double_summation()*s1p1z_z.double_summation();
   cout << "s1p1_p3:  " << s1p1z_x.double_summation()*s1p1z_y.double_summation()*s1p1z_z.double_summation();
   cout << "\n\n";


    // Analytical Integral 7.60173450533140338e-01 for s1s2

   //x
//    float L_A_s1_x=0;
   float L_B_s2_x=0;
//    float R_A_s1_x=0;
   float R_B_s2_x=1;
//    float alpha_s1_x=1;
   float beta_s2_x=1;
   S_AB s1s2_x(L_A_s1_x,L_B_s2_x,R_A_s1_x,R_B_s2_x,alpha_s1_x,beta_s2_x);

   //x tot
   cout << "s1s2  x  " << s1s2_x.double_summation();
   cout << "\n\n";
   
   //y
   float L_A_s1_y=0;
   float L_B_s2_y=0;
   float R_A_s1_y=0;
   float R_B_s2_y=1;
   float alpha_s1_y=1;
   float beta_s2_y=1;
   S_AB s1s2_y(L_A_s1_y,L_B_s2_y,R_A_s1_y,R_B_s2_y,alpha_s1_y,beta_s2_y);

   //y tot
   cout << "s1s2 y  " << s1s2_y.double_summation();
   cout << "\n\n";
   
   
   //z
   float L_A_s1_z=0;
   float L_B_s2_z=0;
   float R_A_s1_z=0;
   float R_B_s2_z=2;
   float alpha_s1_z=1;
   float beta_s2_z=1;
   S_AB s1s2_z(L_A_s1_z,L_B_s2_z,R_A_s1_z,R_B_s2_z,alpha_s1_z,beta_s2_z);

   //z tot
   cout << "s1s2 z  " << s1s2_z.double_summation();
   cout << "\n\n";

   //xyz tot
   float s1s2_xyz_tot = s1s2_z.double_summation()*s1s2_y.double_summation()*s1s2_x.double_summation();
   cout << "s1s2 tot xyz: " << s1s2_xyz_tot;
   cout << "\n\n";


    // Analytical Integral -3.80086725266570169e-01 for s1p2 

    //x
   float L_Ap2x_x=1;
//    float L_Bs_x=0;
   float R_Ap2x_x=1;
//    float R_Bs_x=0;
   float alpha_p2x_x=1;
//    float beta_s_x=1;
   S_AB s1p2x_x(L_Ap2x_x,L_Bs_x,R_Ap2x_x,R_Bs_x,alpha_p2x_x,beta_s_x);
   cout << "x "<< s1p2x_x.double_summation();
   cout << "\n\n";
   
   
   float L_Ap2x_y=0;
//    float L_Bs_y=0;
   float R_Ap2x_y=1;
//    float R_Bs_y=0;
   float alpha_p2x_y=1;
//    float beta_s_y=1;
   S_AB s1p2x_y(L_Ap2x_y,L_Bs_y,R_Ap2x_y,R_Bs_y,alpha_p2x_y,beta_s_y);
   cout << "y "<< s1p2x_y.double_summation();
   cout << "\n\n";
   
   
   float L_Ap2x_z=0;
//    float L_Bs_z=0;
   float R_Ap2x_z=2;
//    float R_Bs_z=0;
   float alpha_p2x_z=1;
//    float beta_s_z=1;
   S_AB s1p2x_z(L_Ap2x_z,L_Bs_z,R_Ap2x_z,R_Bs_z,alpha_p2x_z,beta_s_z);
   cout << "z "<< s1p2x_z.double_summation();
   cout << "\n\n";
   
   //x tot
   float s1p2_x_tot=s1p2x_x.double_summation()*s1p2x_y.double_summation()*s1p2x_z.double_summation();
   cout << "s1p2_p1 :  " << s1p2x_x.double_summation()*s1p2x_y.double_summation()*s1p2x_z.double_summation();
   cout << "\n\n";

   //y
   float L_Ap2y_x=0;
//    float L_Bs_x=0;
   float R_Ap2y_x=1;
//    float R_Bs_x=0;
   float alpha_p2y_x=1;
//    float beta_s_x=1;
   S_AB s1p2y_x(L_Ap2y_x,L_Bs_x,R_Ap2y_x,R_Bs_x,alpha_p2y_x,beta_s_x);
   cout << "x "<< s1p2y_x.double_summation();
   cout << "\n\n";
   
   
   float L_Ap2y_y=1;
//    float L_Bs_y=0;
   float R_Ap2y_y=1;
//    float R_Bs_y=0;
   float alpha_p2y_y=1;
//    float beta_s_y=1;
   S_AB s1p2y_y(L_Ap2y_y,L_Bs_y,R_Ap2y_y,R_Bs_y,alpha_p2y_y,beta_s_y);
   cout << "y "<< s1p2y_y.double_summation();
   cout << "\n\n";
   
   
   float L_Ap2y_z=0;
//    float L_Bs_z=0;
   float R_Ap2y_z=2;
//    float R_Bs_z=0;
   float alpha_p2y_z=2;
//    float beta_s_z=1;
   S_AB s1p2y_z(L_Ap2y_z,L_Bs_z,R_Ap2y_z,R_Bs_z,alpha_p2y_z,beta_s_z);
   cout << "z "<< s1p2y_z.double_summation();
   cout << "\n\n";
   
   //y tot
   float s1p2_y_tot=s1p2y_x.double_summation()*s1p2y_y.double_summation()*s1p2y_z.double_summation();
   cout << "s1p2_p2 :  " << s1p2y_x.double_summation()*s1p2y_y.double_summation()*s1p2y_z.double_summation();
   cout << "\n\n";


   //z
   float L_Ap2z_x=0;
//    float L_Bs_x=0;
   float R_Ap2z_x=1;
//    float R_Bs_x=0;
   float alpha_p2z_x=1;
//    float beta_s_x=1;
   S_AB s1p2z_x(L_Ap2z_x,L_Bs_x,R_Ap2z_x,R_Bs_x,alpha_p2z_x,beta_s_x);
   cout << "x "<< s1p2z_x.double_summation();
   cout << "\n\n";
   
   
   float L_Ap2z_y=0;
//    float L_Bs_y=0;
   float R_Ap2z_y=1;
//    float R_Bs_y=0;
   float alpha_p2z_y=1;
//    float beta_s_y=1;
   S_AB s1p2z_y(L_Ap2z_y,L_Bs_y,R_Ap2z_y,R_Bs_y,alpha_p2z_y,beta_s_y);
   cout << "y "<< s1p2z_y.double_summation();
   cout << "\n\n";
   
   
   float L_Ap2z_z=1;
//    float L_Bs_z=0;
   float R_Ap2z_z=2;
//    float R_Bs_z=0;
   float alpha_p2z_z=1;
//    float beta_s_z=1;
   S_AB s1p2z_z(L_Ap2z_z,L_Bs_z,R_Ap2z_z,R_Bs_z,alpha_p2z_z,beta_s_z);
   cout << "z "<< s1p2z_z.double_summation();
   cout << "\n\n";
   
   //z tot
   float s1p2_z_tot=s1p2z_x.double_summation()*s1p2z_y.double_summation()*s1p2z_z.double_summation();
   cout << "s1p3 :  " << s1p2z_x.double_summation()*s1p2z_y.double_summation()*s1p2z_z.double_summation();
   cout << "\n\n";




   //p-p shell
   //p1p1 
   S_AB p1p1_x(1,1,0,0,1,1);
   float p1p1_x_=p1p1_x.double_summation();
   S_AB p1p1_y(0,0,0,0,1,1);
   float p1p1_y_=p1p1_y.double_summation();
   S_AB p1p1_z(0,0,0,0,1,1);
   float p1p1_z_=p1p1_z.double_summation();
   float p1p1_xyz= p1p1_x_*p1p1_y_*p1p1_z_;

   //p1p2
   S_AB p1p2_x(1,0,0,0,1,1);
   float p1p2_x_=p1p2_x.double_summation();
   S_AB p1p2_y(0,1,0,0,1,1);
   float p1p2_y_=p1p2_y.double_summation();
   S_AB p1p2_z(0,0,0,0,1,1);
   float p1p2_z_=p1p2_z.double_summation();
   float p1p2_xyz= p1p2_x_*p1p2_y_*p1p2_z_;

   //p1p3
   S_AB p1p3_x(1,0,0,0,1,1);
   float p1p3_x_=p1p3_x.double_summation();
   S_AB p1p3_y(0,0,0,0,1,1);
   float p1p3_y_=p1p3_y.double_summation();
   S_AB p1p3_z(0,1,0,0,1,1);
   float p1p3_z_=p1p3_z.double_summation();
   float p1p3_xyz= p1p3_x_*p1p3_y_*p1p3_z_;

   //p2p1
   S_AB p2p1_x(0,1,0,0,1,1);
   float p2p1_x_=p2p1_x.double_summation();
   S_AB p2p1_y(1,0,0,0,1,1);
   float p2p1_y_=p2p1_y.double_summation();
   S_AB p2p1_z(0,0,0,0,1,1);
   float p2p1_z_=p2p1_z.double_summation();
   float p2p1_xyz= p2p1_x_*p2p1_y_*p2p1_z_;

   //p2p2
   S_AB p2p2_x(0,0,0,0,1,1);
   float p2p2_x_=p2p2_x.double_summation();
   S_AB p2p2_y(1,1,0,0,1,1);
   float p2p2_y_=p2p2_y.double_summation();
   S_AB p2p2_z(0,0,0,0,1,1);
   float p2p2_z_=p2p2_z.double_summation();
   float p2p2_xyz= p2p2_x_*p2p2_y_*p2p2_z_;
   
   //p2p3
   S_AB p2p3_x(0,0,0,0,1,1);
   float p2p3_x_=p2p3_x.double_summation();
   S_AB p2p3_y(1,0,0,0,1,1);
   float p2p3_y_=p2p3_y.double_summation(); 
   S_AB p2p3_z(0,1,0,0,1,1);
   float p2p3_z_=p2p3_z.double_summation();
   float p2p3_xyz= p2p3_x_*p2p3_y_*p2p3_z_;


   //p3p1
   S_AB p3p1_x(0,1,0,0,1,1);
   float p3p1_x_=p3p1_x.double_summation();
   S_AB p3p1_y(0,0,0,0,1,1);
   float p3p1_y_=p3p1_y.double_summation();
   S_AB p3p1_z(1,0,0,0,1,1);
   float p3p1_z_=p3p1_z.double_summation();
   float p3p1_xyz= p3p1_x_*p3p1_y_*p3p1_z_;

   //p3p2
   S_AB p3p2_x(0,0,0,0,1,1);
   float p3p2_x_=p3p2_x.double_summation();
   cout << " p3p2_x_ " << p3p2_x_;
   cout << "\n\n";
   S_AB p3p2_y(0,1,0,0,1,1);
   float p3p2_y_=p3p2_y.double_summation();
   cout << " p3p2_y_ " << p3p2_y_;
   cout << "\n\n";
   S_AB p3p2_z(1,0,0,0,1,1);
   float p3p2_z_=p3p2_z.double_summation();
   cout << " p3p2_z_ " << p3p2_z_;
   cout << "\n\n";
   float p3p2_xyz= p3p2_x_*p3p2_y_*p3p2_z_;
   
   //p3p3
   S_AB p3p3_x(0,0,0,0,1,1);
   float p3p3_x_=p3p3_x.double_summation();
   cout << " p3p3_x_ " << p3p3_x_;
   cout << "\n\n";
   S_AB p3p3_y(0,0,0,0,1,1);
   float p3p3_y_=p3p3_y.double_summation();
   cout << " p3p3_y_ " << p3p3_y_;
   cout << "\n\n";
   S_AB p3p3_z(1,1,0,0,1,1);
   float p3p3_z_=p3p3_z.double_summation();
   cout << " p3p3_z_ " << p3p3_z_;
   cout << "\n\n";

   float p3p3_xyz= p3p3_x_*p3p3_y_*p3p3_z_;


   arma::mat init = {{p1p1_xyz, p1p2_xyz,p1p3_xyz},{p2p1_xyz,p2p2_xyz,p2p3_xyz},{p3p1_xyz,p3p2_xyz,p3p3_xyz}};
   init.print("init matrix");





   /*Assemble the set of integrals for a shell pair by calling your S ABx routine to evaluate S ABy and S ABz and forming the final product shown in Eq. 2.3. There will be 1, 3, 9 numbers for s −s, p −s, p −p shell pairs. Compare with what your class-mates are getting and what Jiashu provides you with in lab! */
   return 0;
}