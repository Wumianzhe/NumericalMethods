#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define avg(a, b) (((a) + (b)) / 2)
#define eps_low 1
#define eps_high 14

typedef struct result {
  double x;
  int iter_count;
} result_t;

result_t secant(double(f)(double), double x_0, double x_1, double m_1,
                double m_2, double eps);
result_t bisection(double(f)(double), double a, double b, double eps);
double poly(double x) {
  // 2x^4+8x^3+8x^2-1
  return 2 * x * x * x * x + 8 * x * x * x + 8 * x * x - 1;
  /* return ((2 * x + 8) * x + 8) * x * x -1; */
}
double transc(double x) {
  // (x-3)cos(x)=1
  return (x - 3) * cos(x) - 1;
}

int main(int argc, char *argv[]) {
  const int count = eps_high - eps_low + 1;
  result_t results[2][count];
  FILE *file = fopen("results.csv", "w");
  char *format;
  if (file == stdout) {
    format = "e: %e, result: %0.16lf, f(x): %0.12le, iterations: %d\n";
  } else {
    format = "%e, %0.16lf, %0.12le, %d\n";
  }
  for (int i = 0; i < count; i++) {
    results[0][i] = bisection(poly, 0, 1, pow(10, -eps_low - i));
    results[1][i] = bisection(transc, 5, 5.4, pow(10, -eps_low - i));
  }
  /* fprintf(file,"Polynominal\n"); */
  for (int i = 0; i < count; i++) {
    fprintf(file, format, pow(10, -eps_low - i), results[0][i].x,
            poly(results[0][i].x), results[0][i].iter_count);
  }
  /* fprintf(file,"Transcendental\n"); */
  for (int i = 0; i < count; i++) {
    fprintf(file, format, pow(10, -eps_low - i), results[1][i].x,
            transc(results[1][i].x), results[1][i].iter_count);
  }
  /* fprintf(file,"Secant Method\n"); */
  for (int i = 0; i < count; i++) {
    results[0][i] = secant(poly, 1.84, 1.8, 2, 88, pow(10, -eps_low - i));
    // with 3,4 falls down to -1.7814...
    results[1][i] =
        secant(transc, 5.4, 5.35, 2.2, 1.351, pow(10, -eps_low - i));
  }
  /* fprintf(file,"Polynominal\n"); */
  for (int i = 0; i < count; i++) {
    fprintf(file, format, pow(10, -eps_low - i), results[0][i].x,
            poly(results[0][i].x), results[0][i].iter_count);
  }
  /* fprintf(file,"Transcendental\n"); */
  for (int i = 0; i < count; i++) {
    fprintf(file, format, pow(10, -eps_low - i), results[1][i].x,
            transc(results[1][i].x), results[1][i].iter_count);
  }
  fclose(file);
  // execute octave as subprocess for fzero data
  execl("/usr/bin/octave", "octave", "lab_fzero.m", (char *)NULL);
  return 0;
}

result_t bisection(double(f)(double), double a, double b, double eps) {
  result_t res = {0, 0};
  while (fabs(a - b) > 2 * eps) {
    double c = avg(a, b);
    if (f(a) * f(c) < 0) {
      b = c;
    } else {
      a = c;
    }
    res.iter_count++;
  }
  res.x = avg(a, b);
  return res;
}

result_t secant(double(f)(double), double x_0, double x_1, double m_1,
                double m_2, double eps) {
  result_t res = {0, 0};
  double prev, cur = x_0, next = x_1;
  // m_1 --- min of first derivative
  // m_2 --- max of second derivative

  // for numbers far away from 0 just eps might be less than machine epsilon
  do {
    prev = cur;
    cur = next;
    // can lead to division by 0 in later iterations (iteration 11)
    next = cur + (cur - prev) / (f(prev) / f(cur) - 1);
    // both can do so
    /* next = cur+f(cur)*((cur-prev)/(f(prev)-f(cur))); */
    res.iter_count++;
  } while (fabs(pow(fabs(next - cur), 1.61) * m_2 / 2 / m_1) > eps);
  // while (fabs(next-cur)>fabs(eps*next));
  res.x = next;
  return res;
}
