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
    cout << "\n\n";
    cout << "The unnormaizlied matrix for H2";
    cout << "\n\n";
    cout << H2;
    cout << "\n\n";
     
    arma::mat S_12= orth_trans(H2);
    cout << "The S^-1/2 matrix ";
    cout << "\n\n";
    cout << S_12; 
    cout << "\n\n";

    arma::mat diag_h = {{-13.6,0}, {0,-13.6}};
    cout << "\n\n";
    cout << "The H matrix with assoc ionization potential";
    cout << "\n\n";
    cout << diag_h; 
    cout << "\n\n";

    arma::mat hmat_get = h_mat_get(diag_h,H2);
    cout << "\n\n";
    cout << "The H matrix with the offdiagonal calculated";
    cout << "\n\n";
    cout<< hmat_get; 
    cout << "\n\n";

    arma::mat H_mat = {{-13.6,-14.81},{-14.81,-13.6}};
    cout << "\n\n";
    cout << "The actual H matrix with paramateres manually input in";
    cout << "\n\n";
    cout << H_mat; 
    cout << "\n\n";

    arma::mat curly_H=ham_ortho(H_mat,S_12);
    cout << "\n\n";
    cout << "The H prime matrix to input in the next step (diag/orthogonalization) ";
    cout << "\n\n";
    cout<< curly_H; 
    cout << "\n\n";

    arma::mat vprim=h_diag(curly_H,S_12);
    cout << "\n\n";
    cout << "The eigenvalue/eigenvector solver where the eigenvector solved is aka C_mat";
    cout << "\n\n";
    cout <<vprim; 
    cout << "\n\n";

    arma::mat c = get_v(H2);
    cout << "\n\n";
    cout << "The MO Overlap orbital is just the overlap orbital transposed basically (S^-1/2 * S^1/2) ";
    cout << "\n\n"; 
    cout << c;
    cout << "\n\n";

    arma::mat tot_eng = total_energy(curly_H,S_12);
    cout << "\n\n";
    cout << "The epsilon eigenvalues summed up to give total energy of H2 ";
    cout << "\n\n";
    cout << tot_eng; 
    cout << "\n\n";
    
    arma::mat tot_eng_sep = total_energy_sep(curly_H,S_12);
    cout << "\n\n";
    cout << "The epsilon eigenvalues not summed to give antibonding MO energy and bonding MO energy ";
    cout << "\n\n";
    cout << tot_eng_sep; 
    cout << "\n\n";
    
    arma::mat diff = diff_calc(tot_eng_sep,tot_eng);
    cout << "\n\n";
    cout << "The epsilon eigenvalue for the ground state H2 MO bonding subtracted by the total energy  ";
    cout << "\n\n";
    cout << diff; 
    cout << "\n\n";

}


