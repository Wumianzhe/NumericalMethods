#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double a = -2.0;
const double b = -0.3;

double f(double x) { return x * x * x * x * x - 2.2 * x * x * x + 0.5 * x * x - 7 * x - 3.4; }
// double f(double x) {return cos(x);}
double F(double x) { return pow(x, 6) / 6 - 2.2 / 4 * pow(x, 4) + 0.5 / 3 * pow(x, 3) - 3.5 * x * x - 3.4 * x; }

double integrate(double a, double b, double (*f)(double), int N);

int main(int argc, char* argv[]) {
    double epsLow = 1;
    double epsHigh = 8;
    double precise = F(b) - F(a);
    for (double eps = pow(10, -epsLow); eps >= pow(10, -epsHigh); eps /= 10) {
        int k = 1;
        double prev = integrate(a, b, f, k);
        double next = integrate(a, b, f, ++k);
        while (fabs(next - prev) / (2 * 2 - 1) > eps) {
            prev = next;
            next = integrate(a, b, f, ++k);
        }
        /* file << k << ", " << abs(precise - next) << ", " << eps << "\n"; */
        printf("%e, %d, %e\n", eps, k, fabs(precise - next));
    }

    return 0;
}

double integrate(double a, double b, double (*f)(double), int k) {
    int N = pow(2, k - 1) + 1;
    double h = (b - a) / N;
    double sum = -f(a) / 2 - f(b) / 2;
    for (int i = 0; i < N; i++) {
        sum += f(a + h * i);
    }
    sum *= h;
    return sum;
}
