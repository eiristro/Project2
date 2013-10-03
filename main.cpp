#include <iostream>
#include <cmath>
#include "armadillo"
#include <lib.cpp>

using namespace std;
using namespace arma;

double maxoffdiag(mat A, int& k, int& l, int n) {
    // Algorithm that finds the element above the diagonal
    // with the largest absolute value and returns the
    // coordinates (k, l)
    double max_offdiag = 0.0;

    for (int i = 0; i<n; i++){
        for (int j = i + 1; j<n; j++) {
            if (fabs(A(i, j)) > max_offdiag) {
                max_offdiag = fabs(A(i, j));
                k = i;
                l = j;
            }
        }
    }
    return max_offdiag;
}

void Jacobi_rotation(mat& A, int k, int l, int n) {
    // Algorithm running one Jacobi rotation on armadillo matrix A
    // and changing it to the similar matrix B = S_dagger * A * S

    // First we find the values of t, tau, s and c
    double s, c, t;
    if (A(k, l) != 0) {
        double tau;
        tau = (A(l, l) - A(k, k))/(2.0*A(k, l));
        if (tau > 0) {
            t = -tau + sqrt(1 + tau*tau);
        }
        else {
            t = -tau - sqrt(1 + tau*tau);
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    double a_ll, a_kk, a_ik, a_il;

    // Find the elements with index k and l
    a_ll = A(l, l);
    a_kk = A(k, k);
    A(k, k) = c*c*a_kk - 2.0*c*s*A(k, l) + s*s*a_ll;
    A(l, l) = c*c*a_ll + 2.0*c*s*A(k, l) + s*s*a_kk;
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    // Then the rest of the elements
    for (int i = 0; i<n; i++) {
        if (i != k && i != l) {
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = c*a_ik - s*a_il;
            A(k, i) = A(i, k);
            A(i, l) = c*a_il + s*a_ik;
            A(l, i) = A(i, l);
        }
    }
}


int main()
{
    int n, n_step;
    double h, rho_min, rho_max;

    // Defining the constants. rho_max should ideally be infinite
    // but for the rotation to find the eigenvalues it has to be
    // considerably lower than n. n is the size of A while n_step
    // is the number of steps from rho_min to rho_max
    rho_min = 0;
    rho_max = 5;
    n = 200;
    n_step = n + 1;
    h = (rho_max - rho_min)/n_step;

    // Setting up the potential V
    colvec V(n_step);
    for (int i = 0; i < n_step; i++) {
        V(i) = pow((rho_min + i*h), 2.0);
    }

    // d and e are the diagonal and off-diagonal elements of a
    // used in the function tqli from lib.cpp
    double d[n], e[n];

    // Setting up the matrix running from V(1) to V(n_step - 1)
    mat A(n, n);
    for (int i = 1; i<n_step; i++) {
        for (int j = 1; j < n_step; j++) {
            if (i == j) {
                A(i-1, j-1) = 2.0/(h*h) + V(i);
                d[i-1] = A(i-1, j-1); }
            else if ((j == i+1) || (j == i-1)) {
                A(i-1, j-1) = -1.0/(h*h);
                e[i-1] = A(i-1, j-1); }
            else {
                A(i-1, j-1) = 0;
            }
        }
    }

    // Initializing the variables controlling the while loop
    int k, l;
    double epsilon = 1.0e-8;
    double max_num_iter = (double) n * (double) n * (double) n;
    int iter = 0;
    double max_offdiag = maxoffdiag(A, k, l, n);

    while (fabs(max_offdiag) > epsilon && (double) iter < max_num_iter) {
        max_offdiag = maxoffdiag(A, k, l, n);
        Jacobi_rotation(A, k, l, n);
        iter++;
    }

    // Making a vector with the diagonal (eigenvalue) elements
    // this can then be sorted by the armadillo function sort to
    // give the eigenvaules by size
    colvec a(n);
    for (int i = 0; i<n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                a(i) = A(i, i); }
            }
        }

    a = sort(a);

    double **z;
    z = new double * [n];
    for (int i = 0; i < n; i++) {
        z[i] = new double[n];
    }
    tqli(d, e, n, z);

    cout << "Number of iterations: " << iter << endl;
    cout << "rho_ max: " << rho_max << ", n: " << n << "\n" << endl;
    for (int i = 0; i<5; i++) {
        cout << d[i] << ", " << a(i) << endl;
    }

    return 0;
}
