#include "matrix.h"
#include <cmath>
#include <fstream>
#include <functional>
#include <ios>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;
struct mesh_t {
    mesh_t(int size);
    inline double h(int i) const { return x[i] - x[i - 1]; }
    vector<double> x;
    vector<double> y;
    int size;
    mesh_t operator=(const vector<pair<double, double>>& m) {
        size = m.size();
        x.resize(size);
        y.resize(size);
        for (int i = 0; i < size; i++) {
            x[i] = m[i].first;
            y[i] = m[i].second;
        }
        return *this;
    }
};
double newtonApprox(mesh_t& mesh, vector<double>& diffs, double x);
double maxDiv(const mesh_t& mesh, function<double(double)> f);
vector<double> divDiff(mesh_t vals);
mesh_t meshGen(double a, double b, int n, function<double(double)> f);
matrix_t coefMat(const mesh_t& mesh, function<double(double)> df);
matrix_t constsMat(const matrix_t& M, const mesh_t& mesh);
double splineApprox(const matrix_t& CM, const double x, const mesh_t& mesh);
matrix_t ThomasRL(const matrix_t B);

const double a = 1;
const double b = 2.5;
const double c = 1.7;

const vector<int> meshSizes = {3, 5, 7};
const vector<int> kVals = {0, 3, 7};
const int graphSize = 1000;
const int nlow = 3;
const int nhigh = 50;

double f(double x) { return 1.0 / tan(x) + x * x; }
double df(double x) { return -1.0 / (sin(x) * sin(x)) + 2 * x; }
// double df(double x) { return 0; }

double g(double x) { return f(c) + abs(f(x) - f(c)); }
double dg(double x) { return df(x) * ((x < c) ? (-1) : (1)); }

void splineSuite(const string& dir, const string& name, function<double(double)> f, function<double(double)> df);
void newtonSuite(const string& dir, const string& name, function<double(double)> f);

int main(int argc, char* argv[]) {
    string dir = argv[0];
    dir = dir.substr(0, dir.find_last_of('/') + 1);

    splineSuite(dir, "fspline", f, df);
    splineSuite(dir, "gspline", g, dg);
    newtonSuite(dir, "fpoly", f);
    newtonSuite(dir, "gpoly", g);
    return 0;
}

vector<double> divDiff(mesh_t vals) {
    int size = vals.size;
    matrix_t coef(size, size);
    for (int i = 0; i < size; i++) {
        coef(i, 0) = vals.y[i];
    }
    for (int i = 1; i < size; i++) {
        for (int j = 0; j < size - i; j++) {
            coef(j, i) = (coef(j, i - 1) - coef(j + 1, i - 1)) / (vals.x[j] - vals.x[j + i]);
        }
    }
    vector<double> ret(size);
    for (int i = 0; i < size; i++) {
        ret[i] = coef(0, i);
    }
    return ret;
}

double newtonApprox(mesh_t& mesh, vector<double>& diffs, double x) {
    int size = mesh.size;
    double res = diffs[0];
    for (int i = 1; i < size; i++) {
        double prod = 1;
        for (int j = 0; j < i; j++) {
            prod *= x - mesh.x[j];
        }
        res += diffs[i] * prod;
    }
    return res;
}

double maxDiv(mesh_t& mesh, function<double(double)> f) {
    int size = mesh.size;
    double d = 0;
    auto dd = divDiff(mesh);
    for (int i = 0; i < size - 1; i++) {
        double avg = (mesh.x[i] + mesh.x[i + 1]) / 2;
        double pval = newtonApprox(mesh, dd, avg);
        d = max(d, abs(pval - f(avg)));
    }
    return d;
}

mesh_t::mesh_t(int _size) : size(_size) {
    x.resize(_size);
    y.resize(_size);
}
matrix_t coefMat(const mesh_t& mesh, function<double(double)> df) {
    int size = mesh.size;
    matrix_t mat(size, 4);
    // -2*h_1*M_0 - h_1*M_1 = 6(f'(x_0) - (y_1-y_0)/h_1)
    mat(0, 0) = 0;
    mat(0, 1) = -(2.0 / 6) * mesh.h(1);
    mat(0, 2) = -(1.0 / 6) * mesh.h(1);
    mat(0, 3) = df(mesh.x[0]) - (mesh.y[1] - mesh.y[0]) / mesh.h(1);
    for (int i = 1; i < size - 1; i++) {
        double h = mesh.h(i);
        double h_p = mesh.h(i + 1);
        mat(i, 0) = h / (h + h_p);
        mat(i, 1) = 2;
        mat(i, 2) = h_p / (h + h_p);
        mat(i, 3) = 6 / (h + h_p) * ((mesh.y[i + 1] - mesh.y[i]) / h_p - (mesh.y[i] - mesh.y[i - 1]) / h);
    }
    // -4h_n*M_n + h_n M_{n-1} =6*(f'(x_n)-(y_n-y_{n-1})/h_n)
    int n = size - 1;
    mat(n, 1) = (2.0 / 6) * mesh.h(n);
    mat(n, 2) = (1.0 / 6) * mesh.h(n);
    mat(n, 3) = df(mesh.x[n]) - (mesh.y[n] - mesh.y[n - 1]) / mesh.h(n);
    return mat;
}

mesh_t meshGen(double a, double b, int n, function<double(double)> f) {
    mesh_t mesh(n + 1);
    for (int i = 0; i <= n; i++) {
        mesh.x[i] = a + (b - a) / n * i;
        mesh.y[i] = f(mesh.x[i]);
    }
    return mesh;
}

double splineApprox(const matrix_t& CM, const double x, const mesh_t& mesh) {
    int i = distance(mesh.x.begin(), find_if(mesh.x.begin(), mesh.x.end(), [x](double x_i) { return x_i > x; }));
    if (i == mesh.size) {
        i--;
    }
    double h = mesh.h(i);
    return CM(i - 1, 0) * pow(mesh.x[i] - x, 3) / (6 * h) + CM(i, 0) * pow(x - mesh.x[i - 1], 3) / (6 * h) +
           CM(i, 2) * (x - mesh.x[i - 1]) + CM(i, 1);
}

matrix_t constsMat(const matrix_t& M, const mesh_t& mesh) {
    int size = M.rows;
    matrix_t res(size, 3);
    res(0, 0) = M[0];
    for (int i = 1; i < size; i++) {
        double h = mesh.h(i);
        res(i, 0) = M[i];
        res(i, 1) = mesh.y[i - 1] - M[i - 1] * h * h / 6;
        res(i, 2) = (mesh.y[i] - mesh.y[i - 1]) / h - h / 6 * (M[i] - M[i - 1]);
    }
    return res;
}

void splineSuite(const string& dir, const string& name, function<double(double)> f, function<double(double)> df) {
    matrix_t polyvals(graphSize + 1, meshSizes.size() + 2);
    mesh_t graphMesh = meshGen(a, b, graphSize, f);
    for (int i = 0; i <= graphSize; i++) {
        polyvals(i, 0) = graphMesh.x[i];
        polyvals(i, polyvals.cols - 1) = graphMesh.y[i];
    }
    for (int j = 0; j < meshSizes.size(); j++) {
        mesh_t mesh = meshGen(a, b, meshSizes[j], f);
        matrix_t M = ThomasAlg(coefMat(mesh, df));
        matrix_t CM = constsMat(M, mesh);
        for (int i = 0; i <= graphSize; i++) {
            polyvals(i, j + 1) = splineApprox(CM, graphMesh.x[i], mesh);
        }
    }

    ofstream ofstr(dir + name + ".bin", ios::binary);
    polyvals.write(ofstr);
    ofstr.close();

    matrix_t nvals(nhigh - nlow + 1);
    for (int i = nlow; i <= nhigh; i++) {
        mesh_t graphMesh = meshGen(a, b, i * 10, f);
        mesh_t mesh = meshGen(a, b, i, f);
        matrix_t M = ThomasAlg(coefMat(mesh, df));
        matrix_t CM = constsMat(M, mesh);
        nvals[i - nlow] = accumulate(graphMesh.x.begin(), graphMesh.x.end(), 0.0, [&](const double d, const double x) {
            return max(d, abs(f(x) - splineApprox(CM, x, mesh)));
        });
    }
    ofstr.open(dir + "n" + name + ".bin", ios::binary);
    nvals.write(ofstr);
    ofstr.close();

    // matrix_t dvals((dhigh - dlow) / dstep + 1, 2);
    // mesh_t mesh = meshGen(a, b, 100, f);
    // int i = 0;
    // for (double m = dlow; m <= dhigh + 1e-6; m += dstep) {
    //     matrix_t mat = coefMat(mesh, [=](double x) { return m * df(x); });
    //     matrix_t M = ThomasAlg(mat);
    //     matrix_t CM = constsMat(M, mesh);
    //     dvals(i, 0) = m;
    //     dvals(i++, 1) = accumulate(graphMesh.x.begin(), graphMesh.x.end(), 0.0,
    //                                [&](double d, double x) { return max(d, abs(f(x) - splineApprox(CM, x, mesh)));
    //                                });
    // }
    // ofstr.open(dir + "d" + name + ".bin", ios::binary);
    // dvals.write(ofstr);
    // ofstr.close();
}

void newtonSuite(const string& dir, const string& name, function<double(double)> f) {
    matrix_t polyvals(graphSize + 1, meshSizes.size() + 2);
    mesh_t graphMesh = meshGen(a, b, graphSize, f);
    for (int i = 0; i <= graphSize; i++) {
        polyvals(i, 0) = graphMesh.x[i];
        polyvals(i, polyvals.cols - 1) = graphMesh.y[i];
    }
    for (int j = 0; j < meshSizes.size(); j++) {
        auto mesh = meshGen(a, b, meshSizes[j], f);
        auto dd = divDiff(mesh);
        for (int i = 0; i <= graphSize; i++) {
            polyvals(i, j + 1) = newtonApprox(mesh, dd, graphMesh.x[i]);
        }
    }
    ofstream polystream(dir + name + ".bin", ios::binary);
    polyvals.write(polystream);
    polystream.close();

    matrix_t nvals(nhigh - nlow + 1);
    for (int i = nlow; i <= nhigh; i++) {
        auto mesh = meshGen(a, b, i, f);
        nvals[i - nlow] = maxDiv(mesh, f);
    }
    ofstream nstream(dir + "n" + name + ".bin", ios::binary);
    nvals.write(nstream);
    // ofstream nstreamC(dir + "nfpoly.csv");
    // for (int i = 0; i < nvals.rows; i++) {
    //     for (int j = 0; j < nvals.cols; j++) {
    //         nstreamC << nvals(i, j) << ",";
    //     }
    //     nstreamC << endl;
    // }
    nstream.close();
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
