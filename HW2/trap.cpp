// C++ program to implement Trapezoidal rule
#include<stdio.h>
#include<cmath>
#include<iostream>

using namespace std;
// A sample function whose definite integral's
// approximate value is computed using Trapezoidal
// rule
double y(int x)
{
    // int alpha_=-1;
    // int X_a=1;
    // int X_b=2;
    // int beta_=1;

    return exp(-pow(x-1,2)+pow(x-2,2));
}


// Function to evaluate the value of integral
double trapezoidal(double a, double b, double n)
{
    // Grid spacing
    double h = (b-a)/n;
 
    // Computing sum of first and last terms
    // in above formula
    double s = y(a)+y(b);
 
    // Adding middle terms in above formula
    for (int i = 1; i < n; i++)
        s += 2*y(a+i*h);
 
    // h/2 indicates (b-a)/2n. Multiplying h/2
    // with s.
    return (h/2)*s;
}
 
// Driver program to test above function
int main()
{
    // Range of definite integral
    double x0 = 0;
    double xn = 10;
 
    // Number of grids. Higher value means
    // more accuracy
    int n = 20;
 
    printf("Value of integral is %6.4f\n",
                  trapezoidal(x0, xn, n));
    return 0;
}