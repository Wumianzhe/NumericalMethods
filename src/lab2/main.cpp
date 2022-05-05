#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iostream>
using namespace std;
#include "matrix.h"

inline int sign(double n) { return (n >= 0); }

matrix_t reflectionTransform(matrix_t& B);
matrix_t condGen(double cond, int size);
pair<double, double> reflectionVarCond(double cond, matrix_t X);
matrix_t solveAb(matrix_t& A, matrix_t& b);

int main(int argc, char* argv[]) {
    srand(time(0));

    const int size = 10;
    const double condLow = 1;
    const double condMult = 2;
    const int iter = 30;
    ofstream file("2-cond.csv");

    matrix_t X(size, 1);
    for (int i = 0; i < size; i++) {
        X[i] = mat::random(-10, 10);
    }

    double cond = condLow;
    file << std::scientific;
    for (int i = 0; i < iter; i++) {

        pair<double, double> res = reflectionVarCond(cond, X);
        file << cond << ", " << res.first << ", " << res.second << endl;
        cond *= condMult;
    }
    file.close();
    cout << "Cond done" << endl;

    const int condGood = 10;
    const int condBad = 10000;
    const double distLow = 1e-7;
    const double distHigh = 1e-1;
    const double distStep = 1.5;
    matrix_t AGood = condGen(condGood, size);
    matrix_t ABad = condGen(condBad, size);
    matrix_t bGood = AGood * X;
    matrix_t bBad = ABad * X;
    matrix_t db(size, 1);

    cout << AGood << endl;

    file.open("2-dist.csv");
    file << std::scientific;
    for (double dist = distLow; dist < distHigh; dist *= distStep) {
        for (int i = 0; i < size; i++) {
            db[i] = mat::random(-1, 1);
        }
        // ||db||/||b|| = dist
        matrix_t btGood = bGood + db * (norm(bGood) * dist / norm(db));
        matrix_t btBad = bBad + db * (norm(bBad) * dist / norm(db));

        matrix_t resGood = solveAb(AGood, btGood);
        matrix_t resBad = solveAb(ABad, btBad);
        file << dist << ", " << norm(resGood - X) / norm(X) << ", " << norm(resBad - X) / norm(X) << endl;
    }
    file.close();
    return 0;
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

pair<double, double> reflectionVarCond(double cond, matrix_t X) {
    int size = X.rows;
    matrix_t A = condGen(cond, size);
    // for (int i = 0; i < size; i++) {
    //     X[i] = mat::random(-10, 10);
    // }

    matrix_t b = A * X;
    matrix_t res = solveAb(A, b);
    return {norm(res - X), norm(A * res - b)};
}

matrix_t solveAb(matrix_t& A, matrix_t& b) {
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

matrix_t condGen(double cond, int size) {
    matrix_t mat(size, size);
    // diagonal linspace
    mat(0, 0) = 1;
    for (int i = 1; i < size - 1; i++) {
        mat(i, i) = (cond / size) * (i + 1);
    }
    mat(size - 1, size - 1) = cond;

    // Hausholder matrices
    matrix_t W1(size, 1), W2(size, 1);
    for (int i = 0; i < size; i++) {
        W1[i] = mat::random(0, 1);
        W2[i] = mat::random(0, 1);
    }
    W1 = W1 / W1.norm();
    W2 = W2 / W2.norm();

    matrix_t Q1 = eyes(size) - W1 * W1.T() * 2;
    matrix_t Q2 = eyes(size) - W2 * W2.T() * 2;
    mat = Q1 * mat * Q2;
    mat = mat / sqrt(cond);
    return mat;
}
