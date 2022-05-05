#include "matrix.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>

using namespace std;

inline int sign(double n) { return (n >= 0) ? 1 : -1; }

matrix_t posDefGen(const int size, double& l_min, double& l_max);
double matNorm2(const matrix_t& A);
int optFPI(const matrix_t& A, const matrix_t& b, const double alpha, const double eps, matrix_t& X);

int main(int argc, char* argv[]) {
    srand(time(0));

    const int size = 10;
    const int eps_low = 1;
    const int eps_high = 13;
    double l_min, l_max;
    matrix_t A = posDefGen(size, l_min, l_max);
    double alpha = 2 / (l_min + l_max);

    // cout << A << optAlpha(A) << endl;

    matrix_t X(size, 1);
    for (int i = 0; i < size; i++) {
        X[i] = mat::random(-10, 10);
    }
    matrix_t b = A * X;
    matrix_t X_rec(size);
    ofstream F("3-7.csv");
    F << std::scientific;
    double eps = pow(10, -eps_low);
    for (int i = eps_low; i <= eps_high; i++) {
        for (int i = 0; i < size; i++) {
            X_rec[i] = 0;
        }
        F << eps << ", " << optFPI(A, b, alpha, eps, X_rec) << ", " << norm(X - X_rec) << ", " << norm(A * X_rec - b)
          << endl;
        fprintf(stdout, "%d ", i);
        fflush(stdout);
        eps /= 10;
    }
    printf("\n");
    return 0;
}

double matNorm2(const matrix_t& A) {
    ofstream fout;
    double res;
    fout.open("a.bin", ios::out | ios::binary);
    A.write(fout);
    fout.close();

    pid_t pid = fork();
    if (pid) {
        wait(NULL);
    } else {
        execlp("octave", "octave", "lab3_norm.m", NULL);
    }

    ifstream fin("a.bin", ios::in | ios::binary);
    char* ptr = (char*)&res;
    fin.read(ptr, sizeof(double));
    fin.close();
    return res;
}

matrix_t posDefGen(const int size, double& l_min, double& l_max) {
    matrix_t mat(size, size);
    // double det = 100;
    // double multiplier = pow(det, (2.0 / size) / 2);
    // mat(0, 0) = pow(det, (1.0 / size - 1) / 2);

    for (int i = 0; i < size; i++) {
        // mat(i, i) = mat(i - 1, i - 1) * multiplier;
        mat(i, i) = i + 1;
    }
    l_min = mat(0, 0);
    l_max = mat(size - 1, size - 1);

    // Hausholder matrix
    matrix_t W1(size, 1);
    for (int i = 0; i < size; i++) {
        W1[i] = mat::random(0, 1);
    }
    W1 = W1 / W1.norm();

    matrix_t Q1 = eyes(size) - W1 * W1.T() * 2;
    mat = Q1.T() * mat * Q1;

    return mat;
}

int optFPI(const matrix_t& A, const matrix_t& b, const double alpha, const double eps, matrix_t& X) {
    matrix_t prev = X;

    // B = E
    matrix_t C = eyes(A.rows) - A * alpha;
    matrix_t g = b * alpha;
    double mult = matNorm2(C) / (1 - matNorm2(C));

    int iter = 0;

    do {
        prev = X;
        X = C * prev + g;
        iter++;
    } while (mult * norm(X - prev) > eps);

    return iter;
}
