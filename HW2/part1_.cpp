#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>      

using std::string;
using std::vector;
using namespace std;

struct Quadrature
{
    int n=20;
    virtual float next()=0;
};

template<class T>
struct Trapzd : Quadrature 
{
    double a,b,s;
    T &func;
    Trapzd() {};
    Trapzd(T &funcc, const double aa, const double bb) :
    func(funcc), a(aa), b(bb){}

    float next()
    {   
        float h = (b-a)/n;
        float s = func(a)+func(b);
        for (int i = 1; i < n-1; i++)
        {
            s += 2*func(a+i*h);
        }
        return (h/2)*s;
    }

    float qtrap(const double eps=1.0e-8) 
    {
        const int JMAX=20;
        double s,olds=0.0;
        Trapzd<T> t(func,a,b);

        for (int j=0;j<JMAX;j++) 
        {
            s=t.next();
            if (j>5)
                if (abs(s-olds) < eps*abs(olds)||(s == 0.0 && olds == 0.0)) return s;
            olds=s;
            return olds;
        }
        throw ("Too many steps in routine qtrap");
    };

};

//s1s1
float primitive_gaussian_s1s1(float x){
    float X_b= 0;
    float X_a = 0;
    float alpha_= 1;
    float beta_=1;
    float l_A=0;
    float l_B=0;


    return pow(x-X_a,l_A)*pow(x-X_b,l_B)*exp(-alpha_*pow(x-X_a,2)-beta_*pow(x-X_b,2));

}

//s1p1x
float primitive_gaussian_s1p1x(float x){
    int X_b= 0;
    int X_a = 0;
    int alpha_= 1;
    int beta_=1;
    int l_A=0;
    int l_B=1;
    return pow(x-X_a,l_A)*pow(x-X_b,l_B)*exp(-alpha_*pow(x-X_a,2)-beta_*pow(x-X_b,2));
}



//s1s2
float primitive_gaussian_s1s2(float x){
    int X_b= 0;
    int X_a = 1;
    int alpha_= 1;
    int beta_=1;
    int l_A=0;
    int l_B=0;


    return pow(x-X_a,l_A)*pow(x-X_b,l_B)*exp(-alpha_*pow(x-X_a,2)-beta_*pow(x-X_b,2));

}

//s1p2x 
double primitive_gaussian_s1p2x(float x){
    int X_b= 0;
    int X_a = 1;
    int alpha_= 1;
    int beta_=1;
    int l_A=1;
    int l_B=0;


    return pow(x-X_a,l_A)*pow(x-X_b,l_B)*exp(-alpha_*pow(x-X_a,2)-beta_*pow(x-X_b,2));

}



int main(){
    //s1s1x
    Trapzd s1s1(primitive_gaussian_s1s1,-3,3);
    cout << "s1s1 " << s1s1.qtrap();
    cout << "\n\n";

    //s1p1x
    Trapzd s1p1x(primitive_gaussian_s1p1x,-3,3);
    cout << s1p1x.qtrap();
    cout << "\n\n";


    //s1s2
    Trapzd s1s2(primitive_gaussian_s1s2,-3,3);
    cout << "s1s2 " << s1s2.qtrap();
    cout << "\n\n";

    //s1p2x
    Trapzd s1p2x(primitive_gaussian_s1p2x,-3,3);
    cout << s1p2x.qtrap();
    cout << "\n\n";


}