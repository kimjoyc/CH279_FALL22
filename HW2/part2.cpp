#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>      

using std::string;
using std::vector;
using namespace std;

class S_AB{
private:
	int L_A;
    int L_B;
	int R_A;
	int R_B;
    int alpha_;
    int beta_;

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
        if (n == 0 || n==1)
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
        return sqrt(M_PI/alpha_+beta_);
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


        for(int i = 0; i<=L_A; i++)
        {
            for(int j=0; j<=L_B; j++)
            {

                float fact_diff_A= L_A-i;
                float fact_diff_B= L_B-j;

                float fact_tot_A = factorial(L_A)/factorial(i)*factorial(fact_diff_A);
                float fact_tot_B = factorial(L_A)/factorial(j)*factorial(fact_diff_B);

                float n = i + j - 1; 

                float increment_tot=double_factorial(n)*pow(gaussian_prod_center()-R_A,L_A-i)*pow(gaussian_prod_center()-R_B,L_B-j)/pow(2*(alpha_+beta_),(i+j)/2);
                
                float sum = gaussian_alpha_beta_exp()*exp_prefecator()*fact_diff_A*fact_diff_B*increment_tot;
                total += sum;



            }

        }
        return total;



    }
    

};



int main()
{
   /*Check your analytical integration code for S ABx against numerical integration from the first problem! Check that cases which should give zero do give zero. We will discuss some of these issues in the compute lab.*/
   float L_A=0;
   float L_B=0;
   float R_A=0;
   float R_B=0;
   float alpha_=1;
   float beta_=1;
   S_AB s(L_A,L_B,R_A,R_B,alpha_,beta_);
   cout << s.double_summation();
   cout << "\n\n";

   /*Assemble the set of integrals for a shell pair by calling your S ABx routine to evaluate S ABy and S ABz and forming the final product shown in Eq. 2.3. There will be 1, 3, 9 numbers for s −s, p −s, p −p shell pairs. Compare with what your class-mates are getting and what Jiashu provides you with in lab! */
   return 0;
}