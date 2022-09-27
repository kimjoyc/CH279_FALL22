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
    virtual double next()=0;
};

template<class T>
struct Trapzd : Quadrature 
{
    double a,b,s;
    T &func;
    Trapzd() {};
    Trapzd(T &funcc, const double aa, const double bb) :
    func(funcc), a(aa), b(bb){}

    double next()
    {   
        double h = (b-a)/n;
        double s = func(a)+func(b);
        for (int i = 1; i < n-1; i++)
        {
            s += 2*func(a+i*h);
        }
        return (h/2)*s;
    }

    double qtrap(const double eps=1.0e-10) 
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


int primitive_gaussian(float x){
    int X_b= 0;
    int X_a = 0;
    int alpha_= 1;
    int beta_=1;
    int l_A=0;
    int l_B=0;


    return pow(x-X_a,l_A)*pow(x-X_b,l_B)*exp(-alpha_*pow(x-X_a,2)-beta_*pow(x-X_b,2));

}

int main(){
    Trapzd two_stypes(primitive_gaussian,0,1);
    cout << two_stypes.qtrap();
    cout << "\n\n";

}