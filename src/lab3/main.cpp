#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <vector>

const double a = -2.0;
const double b = -0.3;
using namespace std;

double f(double x) { return x * x * x * x * x - 2.2 * x * x * x + 0.5 * x * x - 7 * x - 3.4; }
// double f(double x) {return cos(x);}
double F(double x) { return pow(x, 6) / 6 - 2.2 / 4 * pow(x, 4) + 0.5 / 3 * pow(x, 3) - 3.5 * x * x - 3.4 * x; }

double integrate(double a, double b, function<double(double)> f, int k);
std::pair<double, int> epsIntegrate(double a, double b, function<double(double)> f, double eps);

int main(int argc, char* argv[]) {
    string dir = argv[0];
    dir = dir.substr(0, dir.find_last_of('/') + 1);

    ofstream file(dir + "results.csv");
    file.precision(8);
    file << scientific;
    double epsLow = 1;
    double epsHigh = 13;
    double precise = F(b) - F(a);
    for (double eps = pow(10, -epsLow); eps >= pow(10, -epsHigh); eps /= 10) {
        auto res = epsIntegrate(a, b, f, eps);
        int k = res.second;
        double result = res.first;
        file << k << ", " << abs(precise - result) << ", " << eps << "\n";
    }
    file.flush();
    file.close();

    ofstream errStep(dir + "steps.csv");
    errStep.precision(8);
    errStep << scientific;
    int steps = 20;
    for (int i = 1; i <= steps; i++) {
        errStep << i << ", " << abs(precise - integrate(a, b, f, i)) << '\n';
    }
    errStep.flush();
    errStep.close();

    return 0;
}

double integrate(double a, double b, function<double(double)> f, int k) {
    int N = pow(2, k - 1);
    double h = (b - a) / N;
    double sum = f(a) / 2 + f(b) / 2;
    // for (double x = a+h; x < b; x+=h)
    for (int i = 1; i < N; i++) {
        sum += f(a + h * i);
    }
    return sum * h;
}

pair<double, int> epsIntegrate(double a, double b, function<double(double)> f, double eps) {
    int k = 1;
    double prev = integrate(a, b, f, k);
    double next = integrate(a, b, f, ++k);
    while (abs(next - prev) / 3 > eps) {
        prev = next;
        next = integrate(a, b, f, ++k);
    }
    return {next, k};
}
