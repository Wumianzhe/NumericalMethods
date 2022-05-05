#include "matrix.h"
#include <functional>
#include <iostream>

using namespace std;

const double a = 1;
const double b = 2;
double p(double x) { return 1 / x; }
double q(double x) { return 2 / x; }
double f(double x) { return 2 * log(x) / x; }
double y(double x) { return log(x); }
double dy(double x) { return 1 / x; }

matrix_t ThomasAppend(const matrix_t B, const double x_n);
matrix_t coefMat(const double a, const double b, const matrix_t& boundCoefs, int n, std::function<double(double)> p,
                 std::function<double(double)> q, std::function<double(double)> f);
matrix_t mirror(matrix_t&& mat);
matrix_t ThomasRL(const matrix_t B);

int main(int argc, char* argv[]) {
    string dir = argv[0];
    dir = dir.substr(0, dir.find_last_of('/') + 1);

    matrix_t bounds(2, 3);
    // a_0y(a) + a_1y'(a) = A
    // y(a) = 0
    bounds(0, 0) = 1;
    bounds(0, 1) = 1;
    bounds(0, 2) = bounds(0, 0) * y(a) + bounds(0, 1) * dy(a);
    // b_0y(b) + b_1y'(b) = B
    // b_1 = 0
    bounds(1, 0) = 1;
    bounds(1, 1) = 0;
    bounds(1, 2) = bounds(1, 0) * y(b) + bounds(1, 1) * dy(b);

    int n = 10;
    vector<int> steps = {1, 2};
    matrix_t buf(n + 1, 4);
    for (int i = 0; i <= n; i++) {
        double x = a + (b - a) / n * i;
        buf(i, 0) = x;
        buf(i, 1) = y(x);
    }
    // double y_b = boundCoefs(1, 2) / boundCoefs(1, 0);
    for (int i = 0; i < steps.size(); i++) {
        matrix_t coefs = coefMat(a, b, bounds, n * steps[i], p, q, f);
        // matrix_t sol = ThomasAppend(coefs, y_b);
        matrix_t sol = ThomasRL(coefs);
        for (int j = 0; j <= n; j++) {
            buf(j, 2 + i) = sol[j * steps[i]];
        }
    }

    ofstream fVals(dir + "vals.bin", ios::binary);
    buf.write(fVals);
    fVals.close();

    int divs = 10;
    ofstream fErrs(dir + "errs.csv");
    fErrs << scientific;
    fErrs.precision(8);
    for (int i = 0; i < divs; i++) {
        matrix_t sol = ThomasRL(coefMat(a, b, bounds, n * pow(2, i), p, q, f));
        matrix_t precise(n * pow(2, i) + 1);
        for (int j = 0; j <= n * pow(2, i); j++) {
            precise[j] = y(a + (b - a) / (n * pow(2, i)) * j);
        }
        fErrs << (b - a) / (n * pow(2, i)) << ", " << infnorm(precise - sol) << endl;
    }
    fErrs.close();
    return 0;
}

matrix_t ThomasAppend(const matrix_t B, const double x_n) {
    matrix_t sol(B.rows + 1);
    matrix_t tmp = ThomasAlg(B);
    for (int i = 0; i < B.rows; i++) {
        sol[i] = tmp[i];
    }
    sol[B.rows] = x_n;
    return sol;
}

matrix_t coefMat(const double a, const double b, const matrix_t& boundCoefs, int n, std::function<double(double)> p,
                 std::function<double(double)> q, std::function<double(double)> f) {
    double h = (b - a) / n;

    double a_0 = boundCoefs(0, 0), a_1 = boundCoefs(0, 1), A = boundCoefs(0, 2);
    double b_0 = boundCoefs(1, 0), B = boundCoefs(1, 2);
    matrix_t M(n + 1, 4);

    auto x = [=](int i) { return a + h * i; };

    double div0 = 2 + h * p(x(1));
    M(0, 0) = 0;
    M(0, 1) = 2 * h * a_0 - a_1 * (3 - (2 - h * p(x(1))) / div0);
    M(0, 2) = (4 - (4 - 2 * h * h * q(x(1))) / div0) * a_1;
    M(0, 3) = 2 * A * h + a_1 * 2 * h * h * f(x(1)) / div0;

    for (int i = 1; i < n; i++) {
        M(i, 0) = 1 - h / 2 * p(x(i));
        M(i, 1) = h * h * q(x(i)) - 2;
        M(i, 2) = 1 + h / 2 * p(x(i));
        M(i, 3) = h * h * f(x(i));
    }
    // b_1 = 0
    M(n, 0) = 0;
    M(n, 1) = 2 * h * b_0;
    M(n, 2) = 0;
    M(n, 3) = 2 * B * h;

    // // почему работает с +=
    // M(n - 1, 3) += M(n - 1, 2) * B / b_0;
    // M(n - 1, 2) = M(n - 1, 1);
    // M(n - 1, 1) = M(n - 1, 0);
    // M(n - 1, 0) = 0;
    return M;
}

matrix_t ThomasRL(const matrix_t B) {
    int size = B.rows;
    matrix_t res(size);
    matrix_t delta(size), lambda(size);
    lambda[size - 1] = B(size - 1, 3) / B(size - 1, 1);
    delta[size - 1] = -B(size - 1, 0) / B(size - 1, 1);
    for (int i = size - 2; i >= 0; i--) {
        double div = B(i, 1) + B(i, 2) * delta[i + 1];
        delta[i] = -B(i, 0) / div;
        lambda[i] = (B(i, 3) - B(i, 2) * lambda[i + 1]) / div;
    }
    res[0] = lambda[0];
    for (int i = 1; i < size; i++) {
        res[i] = delta[i] * res[i - 1] + lambda[i];
    }
    return res;
}
