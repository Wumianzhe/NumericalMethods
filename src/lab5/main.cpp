#include <cmath>
#include <functional>
#include <iostream>
#include <matrix.h>

double f(double x, double y) { return y / x + x * cos(x); }
double sol(double x) { return x * sin(x); }
// double f(double x, double y) { return 3 * x * exp(-x) - (x + 1) / x * y; }
// double sol(double x) { return x * x * exp(-x); }

using namespace std;

double step(double y_p, double x, double h, function<double(double, double)> f);
pair<double, int> epsStep(const double y_p, const double x_p, const double x_n, const double eps,
                          function<double(double, double)> f);
double epsStepR(double y_p, const double x_p, const double x_n, const double eps, function<double(double, double)> f);

int main(int argc, char* argv[]) {
    string dir = argv[0];
    dir = dir.substr(0, dir.find_last_of('/') + 1);

    double a = M_PI_2;
    double b = 2 * M_PI;
    double n = 25;
    double h_0 = (b - a) / n;
    matrix_t buf(n + 1, 4);
    vector<double> steps = {1, 2};
    for (int i = 0; i < buf.rows; i++) {
        buf(i, 0) = a + h_0 * i;
        buf(i, 1) = sol(buf(i, 0));
    }
    buf(0, 2) = sol(a);
    buf(0, 3) = sol(a);
    // for each of step values
    for (int j = 0; j < steps.size(); j++) {
        // for each segment
        double h = h_0 / steps[j];
        for (int i = 0; i < n; i++) {
            double y = buf(i, 2 + j);
            // step through
            for (int k = 0; k < steps[j]; k++) {
                y = step(y, buf(i, 0) + h * k, h, f);
            }
            buf(i + 1, 2 + j) = y;
        }
    }
    ofstream fVals(dir + "vals.bin", ios_base::binary);
    buf.write(fVals);
    fVals.close();

    int epsLow = -1;
    int epsHigh = -15;
    ofstream fErrs(dir + "errs.csv");
    fErrs << scientific;
    fErrs.precision(8);
    matrix_t errs(n + 1, epsLow - epsHigh + 1);
    int k = 0;
    for (double eps = pow(10, epsLow); eps >= pow(10, epsHigh); eps /= 10) {
        double y = buf(0, 1);
        double err = 0;
        int iter = 0;
        errs(0, k) = eps;
        for (int i = 0; i < n; i++) {
            auto res = epsStep(y, buf(i, 0), buf(i + 1, 0), eps, f);
            y = res.first;
            iter = max(iter, res.second);
            // y = epsStepR(y,buf(i,0),buf(1+1,0),eps,f);
            err = max(err, (errs(i + 1, k) = abs(y - buf(i + 1, 1))));
        }
        fErrs << eps << ", " << iter << ", " << err << endl;
        k++;
    }
    ofstream fErrsBin(dir + "errs.bin");
    errs.write(fErrsBin);
    fErrsBin.close();

    ofstream fDelta(dir + "delta.csv");
    fDelta << scientific;
    fDelta.precision(8);
    for (double d = 1e-15; d < 0.1; d *= 1.5) {
        double y = buf(0, 1) * (1 - d);
        double err = 0;
        for (int i = 0; i < n; i++) {
            y = epsStep(y, buf(i, 0), buf(i + 1, 0), 1e-9, f).first;
            err = max(err, abs(y - buf(i + 1, 1)));
        }
        fDelta << d << ", " << err << endl;
    }
    fDelta.close();

    ofstream fStep(dir + "step.csv");
    fStep << scientific;
    fStep.precision(8);
    for (int k = 0; k < 10; k++) {
        double h = h_0 / pow(2, k);
        double y = buf(0, 1);
        double err = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < pow(2, k); j++) {
                y = step(y, buf(i, 0) + h * j, h, f);
            }
            err = max(err, abs(y - buf(i + 1, 1)));
        }
        fStep << h << ", " << err << endl;
    }
    fStep.close();

    // // перепроверка тестового примера
    // double x = a;
    // double y_p = sol(a);
    // double h = 0.2;
    // double k_1 = f(x, y_p);
    // double k_2 = f(x + h / 2, y_p + h * k_1 / 2);
    // double k_3 = f(x + h / 2, y_p + h * k_2 / 2);
    // double k_4 = f(x + h, y_p + h * k_3);
    // cout << step(M_PI_2, M_PI_2, 0.2, f) << endl;
    // cout << k_1 << endl << k_2 << endl << k_3 << endl << k_4 << endl;
    // cout << "h = 0.1" << endl;
    // h = 0.1;
    // k_1 = f(x, y_p);
    // k_2 = f(x + h / 2, y_p + h * k_1 / 2);
    // k_3 = f(x + h / 2, y_p + h * k_2 / 2);
    // k_4 = f(x + h, y_p + h * k_3);
    // cout << k_1 << endl << k_2 << endl << k_3 << endl << k_4 << endl;
    // y_p = y_p + h / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
    // x = x + h;
    // cout << y_p << endl;
    // k_1 = f(x, y_p);
    // k_2 = f(x + h / 2, y_p + h * k_1 / 2);
    // k_3 = f(x + h / 2, y_p + h * k_2 / 2);
    // k_4 = f(x + h, y_p + h * k_3);
    // cout << k_1 << endl << k_2 << endl << k_3 << endl << k_4 << endl;
    // y_p = y_p + h / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
    // x = x + h;
    // cout << y_p << endl;
    return 0;
}

double step(double y_p, double x, double h, function<double(double, double)> f) {
    double k_1 = f(x, y_p);
    double k_2 = f(x + h / 2, y_p + h * k_1 / 2);
    double k_3 = f(x + h / 2, y_p + h * k_2 / 2);
    double k_4 = f(x + h, y_p + h * k_3);
    return y_p + h / 6 * (k_1 + 2 * k_2 + 2 * k_3 + k_4);
}

pair<double, int> epsStep(const double y_p, const double x_p, const double x_n, const double eps,
                          function<double(double, double)> f) {
    double h = x_n - x_p;
    double y = step(y_p, x_p, h, f);
    double y_n = step(y_p, x_p, h / 2, f);
    y_n = step(y_n, x_p + h / 2, h / 2, f);
    int iter = 1;
    while (abs(y - y_n) / 15 >= eps) {
        iter++;
        y = y_n;
        y_n = y_p;
        h /= 2;
        for (int i = 0; i < pow(2, iter); i++) {
            y_n = step(y_n, x_p + h / 2 * i, h / 2, f);
        }
    }
    return {y_n, iter};
}

double epsStepR(double y_p, const double x_p, const double x_n, const double eps, function<double(double, double)> f) {
    double x = x_p;
    double y;
    double y_n;
    double y_tmp;
    do {
        double h = x_n - x;
        y = step(y_p, x, h, f);
        y_tmp = step(y_p, x, h / 2, f);
        y_n = step(y_tmp, x + h / 2, h / 2, f);
        while (abs(y - y_n) / 15 >= eps) {
            y = y_tmp;
            h /= 2;
            y_tmp = step(y_p, x, h / 2, f);
            y_n = step(y_tmp, x + h / 2, h / 2, f);
        }
        y_p = y_n;
        x += h;
    } while (x < x_n);
    return y_n;
}
