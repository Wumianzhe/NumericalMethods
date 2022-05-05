#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>

using namespace std;

inline int sign(double n) { return (n >= 0) ? 1 : -1; }
pair<int, int> optEl(const matrix_t& A);
double sumBotLeft(const matrix_t& A);
matrix_t solveAb(const matrix_t& A, const matrix_t& b);
matrix_t rotTransform(const matrix_t& A, int i, int j);
matrix_t inverseIteration(const matrix_t& A, int& iter, const double l, const double eps);
matrix_t condGen(int size, double cond, matrix_t& eigvalues);
matrix_t reflectionTransform(matrix_t& B);
matrix_t eigVec(const matrix_t& A);

int JacobiMethod(const matrix_t& A, matrix_t& lambda, double eps) {
    matrix_t M = A;
    int iter = 0;
    while (sumBotLeft(M) > eps) {
        pair<int, int> opt = optEl(M);
        M = rotTransform(M, opt.first, opt.second);
        iter++;
    }
    for (int i = 0; i < M.cols; i++) {
        lambda[i] = M(i, i);
    }
    return iter;
}

int main(int argc, char* argv[]) {
    // random repeatability
    srand(time(0));
    int seed = rand() % 1024;
    if (argc > 1) {
        seed = atoi(argv[1]);
    } else {
        cout << "seed: " << seed << endl;
    }
    srand(seed);

    const int size = 10;
    const double cond = 100;
    const int eps_low = 2;
    const int eps_high = 14;

    matrix_t l_star(size), l(size);
    matrix_t A = condGen(size, cond, l_star);
    // l_star.data = {1, 3, 5};
    // for (int i = 0; i < 3; i++) {
    //     A(i, i) = l_star[i];
    // }
    // matrix_t W(size);
    // W.data = {0.8542, 0.4270, 0.2965};
    // W = W / W.norm();
    // matrix_t Q = eyes(3) - W * W.T() * 2;
    // A = Q.T() * A * Q;
    // cout << A;

    matrix_t vec_star = eigVec(A);
    ofstream F("4-7.csv");

    F << std::scientific;
    double eps = pow(10, -eps_low);
    for (int i = eps_low; i <= eps_high; i++) {
        int iter_jac = JacobiMethod(A, l, eps);
        std::sort(l.data.begin(), l.data.end());
        int iter_invit;

        // if eps value is same as for lambda, INVIT may not converge using imprecise eigenvalue
        matrix_t vec = inverseIteration(A, iter_invit, l[0], eps);
        F << eps << ", " << norm(l - l_star) << ", " << iter_jac << ", " << iter_invit << ", ";
        F << std::min(norm(vec - vec_star), norm(vec + vec_star)) << ", " << norm(A * vec - vec * l[0]) << endl;
        eps /= 10;
        cout << i << ' ';
    }
    cout << endl;
    return 0;
}

matrix_t solveAb(const matrix_t& A, const matrix_t& b) {
    int size = A.rows;
    matrix_t B(size, size + 1);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            B(i, j) = A(i, j);
        }
    }
    for (int i = 0; i < size; i++) {
        B(i, size) = b[i];
    }
    return reflectionTransform(B);
}

double sumBotLeft(const matrix_t& A) {
    double sum = 0;
    for (int j = 0; j < A.cols; j++) {
        for (int i = j + 1; i < A.rows; i++) {
            sum += A(i, j) * A(i, j);
        }
    }
    return sum;
}

matrix_t rotTransform(const matrix_t& A, int i, int j) {
    matrix_t B = A;
    double p = 2 * A(i, j);
    double q = A(j, j) - A(i, i);
    double d = sqrt(p * p + q * q);
    double r = abs(q) / (2 * d);
    double c = sqrt(0.5 + r);
    double s = sqrt(0.5 - r) * sign(p * q);

    for (int k = 0; k < A.cols; k++) {
        B(i, k) = c * A(k, i) - s * A(k, j);
        B(j, k) = s * A(k, i) + c * A(k, j);
        if (k == i) {
            B(i, i) = c * c * A(i, i) + s * s * A(j, j) - c * s * p;
            B(i, j) = 0;
        }
        if (k == j) {
            B(j, i) = 0;
            B(j, j) = s * s * A(i, i) + c * c * A(j, j) + c * s * p;
        }
    }
    for (int k = 0; k < A.rows; k++) {
        B(k, i) = B(i, k);
        B(k, j) = B(j, k);
    }
    return B;
}

matrix_t reflectionTransform(matrix_t& B) {
    for (int i = 0; i < B.rows - 1; i++) {
        double s = 0;
        for (int k = i; k < B.rows; k++) {
            s += B(k, i) * B(k, i);
        }

        // w construction
        matrix_t W(B.rows, 1);
        for (int k = 0; k < i; k++) {
            W[k] = 0;
        }
        W(i, 0) = B(i, i) + sign(B(i, i)) * sqrt(s);
        for (int k = i + 1; k < B.rows; k++) {
            W[k] = B(k, i);
        }
        double beta = 1 / (s + abs(B(i, i)) * sqrt(s));

        matrix_t H = eyes(B.rows) - W * W.T() * beta;
        B = H * B;
    }

    // gaussBack
    matrix_t res(B.rows, 1);
    for (int i = B.rows - 1; i >= 0; i--) {
        double s = 0;
        for (int k = i + 1; k < B.rows; k++) {
            s += B(i, k) * res[k];
        }
        double x = (B(i, B.cols - 1) - s) / (B(i, i));
        res[i] = x;
    }
    return res;
}

matrix_t inverseIteration(const matrix_t& A, int& iter, const double l, const double eps) {
    int size = A.rows;
    iter = 0;
    double lmb = l;
    matrix_t vec(size);
    for (int i = 0; i < size; i++) {
        vec[i] = 1;
    }
    while (norm(A * vec - vec * lmb) / norm(vec) > eps) {
        iter++;
        matrix_t vec_old = vec / norm(vec);
        vec = solveAb(A - eyes(size) * lmb, vec_old);
        // уточнение l
        double sum = 0;
        int count = 0;
        for (int i = 0; i < size; i++) {
            if (abs(vec[i]) > eps) {
                count++;
                sum += vec_old[i] / vec[i];
            }
        }
        sum /= count;
        lmb += sum;
    }
    vec = vec / norm(vec);
    return vec;
}

matrix_t condGen(int size, double cond, matrix_t& eigvalues) {
    matrix_t mat(size, size);
    double multiplier = pow(cond, 1.0 / (size - 1));

    mat(0, 0) = pow(cond, (1.0 / size - 1) / 2);
    eigvalues[0] = mat(0, 0);
    for (int i = 1; i < size; i++) {
        mat(i, i) = mat(i - 1, i - 1) * multiplier;
        eigvalues[i] = mat(i, i);
    }

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

pair<int, int> optEl(const matrix_t& A) {
    pair<int, int> res = {0, 0};

    // finding row with maximal summ of non-diagonal elements
    double maxRow = 0;
    for (int i = 0; i < A.rows; i++) {
        double rowSum = 0;
        for (int j = 0; j < A.cols; j++) {
            if (i != j) {
                rowSum += A(i, j) * A(i, j);
            }
        }
        if (rowSum > maxRow) {
            maxRow = rowSum;
            res.first = i;
        }
    }

    // finding maximal non-diagonal element of the row
    double maxAbs = 0;
    int i = res.first;
    for (int j = 0; j < A.cols; j++) {
        if (j != i) {
            if (abs(A(i, j)) > maxAbs) {
                maxAbs = abs(A(i, j));
                res.second = j;
            }
        }
    }

    return res;
}
matrix_t eigVec(const matrix_t& A) {
    fstream file;
    matrix_t res(A.cols);
    char* ptr = (char*)A.data.data();
    file.open("a.bin", ios::out | ios::binary);
    file.write(ptr, sizeof(double) * A.rows * A.cols);
    file.close();

    pid_t pid = fork();
    if (pid) {
        wait(NULL);
    } else {
        execlp("octave", "octave", "lab4_eig.m", NULL);
    }

    ptr = (char*)res.data.data();
    file.open("a.bin", ios::in | ios::binary);
    file.read(ptr, sizeof(double) * A.cols);
    file.close();
    return res;
}
