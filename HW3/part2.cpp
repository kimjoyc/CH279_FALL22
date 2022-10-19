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
//1. make the orthogonalization transformation: X =S−1/2
arma::mat orth_trans(arma::mat Overlap){
    
    arma::vec S_eigval;
    arma::mat S_eigvec;
    eig_sym(S_eigval, S_eigvec, Overlap);  
    return S_eigvec*sqrt(inv(diagmat(S_eigval)))*S_eigvec.t();
}

// X = S^1/2
arma::mat orth_trans_(arma::mat Overlap){
    
    arma::vec S_eigval;
    arma::mat S_eigvec;
    eig_sym(S_eigval, S_eigvec, Overlap);  
    return S_eigvec*sqrt(diagmat(S_eigval))*S_eigvec.t();
}

//ham matrix 
arma::mat h_mat_get(arma::mat h2, arma::mat h){
    arma::vec h_mu_nu= 1.75*0.5*arma::sum(h2.diag(),0);
    h=h_mu_nu[0]*h-h2;
    return h;
}


//. form the hamiltonian in the orthogonalized basis: H=XT HX
arma::mat ham_ortho(arma::mat Hmat,arma::mat S_inv){

    return S_inv.t()*Hmat*S_inv;
}

//3. diagonalize: HV =Vε
arma::mat h_diag(arma::mat h, arma::mat s){
    
    arma::vec epsilon;
    arma::mat v_prim;
    eig_sym(epsilon, v_prim, h);  
    
    arma::mat v = s*v_prim;
    return v;
}


//4. form the MO coefficients: C =XV
// basis vectors form an orthonormal set, the orbital overlap matrix will be the identity matrix. T
arma::mat get_v(arma::mat H){

    return orth_trans(H)*orth_trans_(H);
}

//calc total energy
/*
You evaluate this by
subtracting twice the energy of the H atom (−27.2 eV ) from your computed H2 energy.
This result is also your next debugging case. The result should be around −4 eV (the
minus sign indicates binding).
*/

arma::mat total_energy(arma::mat h, arma::mat s){
    
    arma::vec epsilon;
    arma::mat v_prim;
    eig_sym(epsilon, v_prim, h);    
    return arma::sum(epsilon,0)*2;
}

arma::mat total_energy_sep(arma::mat h, arma::mat s){
    
    arma::vec epsilon;
    arma::mat v_prim;
    eig_sym(epsilon, v_prim, h);    
    return epsilon*2;
}

arma::mat diff_calc(arma::mat h2, arma::mat h_tot){
    return h2[0]-h_tot; 
}


int main() {
    //H2 test case 
    arma::mat H2 = {{1.1522,0.6223},{0.6223,1.1522}};
    arma::mat S_12= orth_trans(H2);
    arma::mat diag_h = {{-13.6,0}, {0,-13.6}};
    arma::mat hmat_get = h_mat_get(diag_h,H2);
    arma::mat H_mat = {{-13.6,-15},{-15,-13.6}};
    arma::mat curly_H=ham_ortho(H_mat,S_12);
    arma::mat vprim=h_diag(curly_H,S_12);
    arma::mat c = get_v(H2);
    arma::mat tot_eng = total_energy(curly_H,S_12);
    arma::mat tot_eng_sep = total_energy_sep(curly_H,S_12);
    arma::mat diff = diff_calc(tot_eng_sep,tot_eng);
    cout << diff;

}


