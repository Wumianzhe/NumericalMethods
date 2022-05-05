#include "matrix.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

using namespace std;
// typedef vector<pair<double, double>> mesh_t;
struct mesh_t {
    mesh_t(int size);
    inline double h(int i) const { return x[i] - x[i - 1]; }
    vector<double> x;
    vector<double> y;
    int size;
};

const double a = 1;
const double b = 2.5;
const double c = 1.7;
double f(double x) { return 1.0 / tan(x) + x * x; }
double df(double x) { return -1.0 / (sin(x) * sin(x)) + 2 * x; }
double g(double x) { return f(c) + abs(f(x) - f(c)); }
double dg(double x) { return df(x) * ((x < c) ? (-1) : (1)); }

mesh_t meshGen(double a, double b, int n, function<double(double)> f);
matrix_t coefMat(const mesh_t& mesh, function<double(double)> df);
matrix_t constsMat(const matrix_t& M, const mesh_t& mesh);
double splineApprox(const matrix_t& CM, const double x, const mesh_t& mesh);
void graphSuite(const string& dir, const string& name, function<double(double)> f, function<double(double)> df);

const vector<int> meshSizes = {3, 7, 10};
const int graphSize = 1000;
const int nlow = 7;
const int nhigh = 100;
const double dlow = -1;
const double dhigh = 2;
const double dstep = 0.001;

int main(int argc, char* argv[]) {
    string dir = argv[0];
    dir = dir.substr(0, dir.find_last_of('\\') + 1);

    graphSuite(dir, "smooth", f, df);
    graphSuite(dir, "corner", g, dg);
    // mesh_t polymesh = meshGen(a, b, 3, f);
    // matrix_t M = ThomasAlg(coefMat(polymesh, df));
    // cout << M << endl;
    // matrix_t CM = constsMat(M, polymesh);
    // cout << CM << endl;

    // Expression x("x");
    // for (int i = 1; i < 4; i++) {
    //     double h = polymesh.h(i);
    //     auto g_i = CM(i - 1, 0) * pow(polymesh.x[i] - x, 3) / (6 * h) +
    //                CM(i, 0) * pow(x - polymesh.x[i - 1], 3) / (6 * h) + CM(i, 2) * (x - polymesh.x[i - 1]) + CM(i,
    //                1);
    //     cout << "g_" << i << ": " << expand(g_i) << endl;
    // }

    return 0;
}

matrix_t coefMat(const mesh_t& mesh, function<double(double)> df) {
    int size = mesh.size;
    matrix_t mat(size, 4);
    // -2*h_1*M_0 - h_1*M_1 = 6(f'(x_0) - (y_1-y_0)/h_1)
    mat(0, 0) = -(2.0 / 6) * mesh.h(1);
    mat(0, 1) = -(1.0 / 6) * mesh.h(1);
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

mesh_t::mesh_t(int _size) : size(_size) {
    x.resize(_size);
    y.resize(_size);
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

void graphSuite(const string& dir, const string& name, function<double(double)> f, function<double(double)> df) {
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

    matrix_t dvals((dhigh - dlow) / dstep + 1, 2);
    mesh_t mesh = meshGen(a, b, 100, f);
    int i = 0;
    for (double m = dlow; m <= dhigh + 1e-6; m += dstep) {
        matrix_t mat = coefMat(mesh, [=](double x) { return m * df(x); });
        matrix_t M = ThomasAlg(mat);
        matrix_t CM = constsMat(M, mesh);
        dvals(i, 0) = m;
        dvals(i++, 1) = accumulate(graphMesh.x.begin(), graphMesh.x.end(), 0.0,
                                   [&](double d, double x) { return max(d, abs(f(x) - splineApprox(CM, x, mesh))); });
    }
    ofstr.open(dir + "d" + name + ".bin", ios::binary);
    dvals.write(ofstr);
    ofstr.close();
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
