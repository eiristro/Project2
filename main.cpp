#include <iostream>
#include <cmath>
#include "armadillo"
#include <lib.cpp>
#include <time.h>

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
    double h, rho_min, rho_max, omega_r, rho;

    // Defining the constants. rho_max should ideally be infinite
    // but for the rotation to find the eigenvalues it has to be
    // considerably lower than n. n is the size of A while n_step
    // is the number of steps from rho_min to rho_max
    rho_min = 0;
    rho_max = 5;
    n = 100;
    n_step = n + 1;
    h = (rho_max - rho_min)/n_step;

    omega_r = 5;

    // Setting up the potential V
    colvec V(n_step), Vr(n_step);
    for (int i = 0; i < n_step; i++) {
        rho = i*h;
        Vr(i) = pow(omega_r*rho, 2) + 1.0/rho;
        V(i) = pow(omega_r*rho, 2);            // Withouth the repulsive force
    }

    // d and e are the diagonal and off-diagonal elements of a
    // used in the function tqli from lib.cpp. The r versions include
    // the repulsive force
    double d[n], dr[n], e[n], er[n];

    // Setting up the matrix running from V(1) to V(n_step - 1)
    mat A(n, n);
    for (int i = 1; i<n_step; i++) {
        for (int j = 1; j < n_step; j++) {
            if (i == j) {
                A(i-1, j-1) = 2.0/(h*h) + V(i);
                d[i-1] = A(i-1, j-1);
                dr[i-1] = 2.0/(h*h) + Vr(i);
            }
            else if ((j == i+1) || (j == i-1)) {
                A(i-1, j-1) = -1.0/(h*h);
                e[i-1] = A(i-1, j-1);
                er[i-1] = e[i-1];
            }
            else {
                A(i-1, j-1) = 0;
            }
        }
    }

    // Finding the time for the Jacobi algo
    clock_t start, finish;
    start = clock();

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

    finish = clock();
    double time_Jacobi = ((double)(finish - start)/CLOCKS_PER_SEC);

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

    // Creating an empty matrix z that will hold the eigenvectors
    // we get from tqli
    double **z;
    z = new double * [n];
    for (int i = 0; i < n; i++) {
        z[i] = new double[n];
    }

    double **zr;
    zr = new double * [n];
    for (int i = 0; i < n; i++) {
        zr[i] = new double[n];
    }

    // Running tqli from lib.cpp which finds the eigenvalues using
    // an alog based on Householders algo, and finding the time it
    // takes to run
    start = clock();

    tqli(d, e, n, z);
    tqli(dr, er, n, zr);

    finish = clock();
    double time_tqli = ((finish - start)/(double)CLOCKS_PER_SEC);

    // Creating an armadillo vector holding the elements of d
    // for easy sorting
    colvec ad(n);
    for (int i = 0; i<n; i++) {
        ad(i) = d[i];
    }
    ad = sort(ad);

    colvec adr(n);
    for (int i = 0; i<n; i++) {
        adr(i) = dr[i];
    }
    adr = sort(adr);

    // Finding the eigenvector corresponding to the lowest eigenvalue
    double min_d = ad(n-1);
    int m;
    for (int i = 0; i<n; i++) {
        if (d[i] < min_d) {
            min_d = d[i];
            m = i;
        }
    }

    colvec min_z(n);
    for (int i = 0; i<n; i++) {
        min_z(i) = z[i][m];
    }


    double min_dr = adr(n-1);
    int mr;
    for (int i = 0; i<n; i++) {
        if (dr[i] < min_dr) {
            min_dr = dr[i];
            mr = i;
        }
    }

    colvec min_zr(n);
    for (int i = 0; i<n; i++) {
        min_zr(i) = zr[i][mr];
    }

    // Outputting relevant variables
    cout << "Number of iterations: " << iter << endl;
    cout << "rho_ max: " << rho_max << ", n: " << n
         << ", omega_r: " << omega_r << endl;
    cout << "Time taken: Jacobi: " << time_Jacobi <<
            ", tqli: " << time_tqli << endl << endl;
    cout << "Non-int,   int,    Jacobi" << endl;
    for (int i = 0; i<5; i++) {
        cout << ad(i) << ", " << adr(i) << ", " << a(i) << endl;
    }

    // Outputting z to a file for plotting in python
    char buffer[50], buffer2[50];
    int omega100 = (int) 100*omega_r;

    cout << omega100;
    sprintf(buffer, "../Project2/zdata%d.dat", omega100);
    min_z.save(buffer, raw_ascii);

    sprintf(buffer2, "../Project2/zdatar%d.dat", omega100);
    min_zr.save(buffer2, raw_ascii);

    return 0;
}
