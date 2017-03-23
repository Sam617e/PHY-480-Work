#include <iostream>
#include <armadillo.h>
#include <iostream>
#include <cmath>
/*
include "jacobi.h"
*/


using namespace arma;
using namespace std;

/*
Jacobi's method for finding eigenvalues
eigenvectors of the symetric matrix A.
The eigenvalues of A will be on the diagonal
of A, with eigenvalue i being A[i][i].
The j-th component of the i-th eigenvector
is stored in R[i][j].
A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n)
n: dimention of matrices
*/

double offdiag(mat &A, int &p,
               int &q, int n){
    double max;
    for (int i = 0; i < n; ++i){
        for ( int j = i+1; j < n; ++j){
            double aij = fabs(A(i,j));
            if ( aij > max){
                max = aij; p = i; q = j;}
        }
    }
return max;}

/*
 This is does the similarity transformation (S^T)A(S) on Matrix A
 converts it into matrix B
*/
mat Jacobi_rotate ( mat &A, mat &R, int k,
                    int l, int n ){
    double s, c;
    if ( A(k,l) != 0.0 ){
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        if ( tau >= 0 ){
            t = 1.0/(tau+ sqrt(1.0 + tau*tau));}
        else{
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));}
        c = 1/sqrt(1+t*t);
        s = c*t;}
    else{
        c = 1.0;
        s = 0.0;}
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for ( int i = 0; i < n; i++ ){
        if ( i != k && i != l ){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);}
    r_ik = R(i,k);
    r_il = R(i,l);
    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;}
return A;
}
// end of function jacobi_rotate

mat test_operation ( mat &A, int n )
{
    for ( int i = 0; i < n; i++)
    {
        //cout << A(i,i) << endl;
        A(i,i) = 1.0;
    }
    //cout << A << endl;
    return A;
}

int main(int argc, char *argv[])
{
    int max = 8;
    vec rho(max, fill::zeros);
    rho(1) = 0;
    rho(max-1) = 7;
    float h = (rho(max-1) - rho(0))/(max);
    for(int i = 1; i < max; i++)
    {
        rho(i) = rho(1) + i*h;
    }

    // Change this value for different strengths
    // of omega_r. We're doing 0.01, 0.5, 1, and 5.
    double freq = 0.01;

    vec V2(max, fill::zeros);
    V2(0) = 0;
    for(int i = 1; i < max; i++)
    {
        V2(i) = rho(i)*rho(i)*freq*freq + 1/rho(i);
    }

    float hh = h*h;

    vec d2(max);
    for(int i = 0; i < max; i++)
    {
        d2(i) = V2(i) + 2/hh;
    }

    float e = -1/hh;

    mat M(max, max, fill::zeros);
    M(max-1,max-1) = d2(max-1);
    for(int i = 0; i < max-1; i++)
    {
        M(i,i) = d2(i);
        M(i,i+1) = e;
        M(i+1,i) = e;
    }
    mat R = eye (max,max);
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, M);

    cout << M << endl;
    cout << R << endl;
    double maxiter = 100; //max * max * max;
    double max_offdiag = 1;//offdiag ( M, p, q, max);

    double tolerance = 1.0E-10;
    int iterations = 0;

    while ( max_offdiag > tolerance &&
            iterations <= maxiter)
    {
    int p, q;
    max_offdiag = offdiag(M, p, q, max);
    cout << max_offdiag << endl << endl;
    Jacobi_rotate(M, R, p, q, max);
    cout << M << endl << endl;
    iterations++;
    }
    cout << M << endl << endl;
    cout << R << endl << endl;
    cout << "Number of iterations: " << iterations << "\n";

    cout << "Armadillo Eigenvalues (Interacting): " << endl << eigval << endl;

    return 0;

}
