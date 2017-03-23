#include <iostream>
#include <armadillo.h>
#include <iostream>
#include <cmath>

using namespace arma;
using namespace std;

/*
Unit Test Program
Will do a test of offdiag for a simple 3x3 matrix
Will also do a test of orthogonality of eigenvectors
*/

double offdiag(mat &A, int &p, int &q, int n){
    double max;
    for (int i = 0; i < n; ++i){
        for ( int j = i+1; j < n; ++j){
            double aij = fabs(A(i,j));
            if ( aij > max){
                max = aij; p = i; q = j;}
        }
    }
return max;}

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


int main(int argc, char *argv[])
{
    mat M(3, 3, fill::zeros);
    for(int i = 0; i < 3; i++){
        M(i,i) = 1;
    }

    M(0,1) = 0;
    M(0,2) = 16;
    M(1,2) = -80;
    M(1,0) = 0;
    M(2,0) = 16;
    M(2,1) = -80;
    cout << "Initial M" << endl;
    cout << M;

    int p;
    int q;
    cout << "Max off-diagonal element is: " << offdiag(M, p, q, 3) << endl;
    mat R = eye (3,3); //identity matrix that stores the original orthonormal eigenvalues
    //The three empty eigenvectors to be computed and dotted later.
    vec v0(3, fill::zeros);
    vec v1(3, fill::zeros);
    vec v2(3, fill::zeros);

    M = Jacobi_rotate(M, R, p, q, 3);
    cout << "M after one rotation: " << endl;
    cout << M;
    cout << "Eigenvectors after one rotation: " << endl;
    cout << R << endl;
    cout << "Test for orthogonality:" << endl;

    for (int i = 0; i < 3; i++){
        v0(i) = R(i,0);
        v1(i) = R(i,1);
        v2(i) = R(i,2);
    }

    cout << "v0 = " << v0;
    cout << "v1 = " << v1;
    cout << "v2 = " << v2;

    cout << "Dot product of v0 and v1 = " << dot(v0, v1) << endl;
    cout << "Dot product of v1 and v2 = " << dot(v1, v2) << endl;
    cout << "Dot product of v2 and v0 = " << dot(v2, v0) << endl << endl;

//Confirming eigenvectors with Armadillo functions:

    mat A(3,3, fill::zeros);

    for (int i = 0; i < 3; i++){
        A(i,0) = v0(i);
        A(i,1) = v1(i);
        A(i,2) = v2(i);
    }
    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, A);
    cout << "Armadillo Eigenvectors:" << endl;
    cout << eigvec;

    return 0;
}
