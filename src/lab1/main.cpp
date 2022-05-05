#include "matrix.h"
#include <cmath>
#include <fstream>
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
double newtonApprox(vector<pair<double, double>> mesh, double x);
double maxDiv(vector<pair<double, double>> mesh, double(f)(double));
vector<double> divDiff(vector<pair<double, double>> vals);
vector<pair<double, double>> meshGen(double a, double b, int n, double(f)(double));

#define c 1.7

double f(double x) { return 1 / tan(x) + x * x; }
double g(double x) { return f(c) + abs(f(x) - f(c)); }
double l(double x) { return 2 * x - log10(x); }

int main(int argc, char* argv[]) {
    const double a = 1;
    const double b = 2.5;
    const vector<int> meshSizes = {3, 5, 7};
    const int graphSize = 100;
    const int nlow = 3;
    const int nhigh = 50;

    string dir = argv[0];
    dir = dir.substr(0, dir.find_last_of('/') + 1);

    // гладкая функция
    matrix_t polyvals(graphSize + 1, meshSizes.size() + 2);
    mesh_t graphMesh(1);
    graphMesh = meshGen(a, b, graphSize, f);
    for (int i = 0; i <= graphSize; i++) {
        polyvals(i, 0) = graphMesh.x[i];
        polyvals(i, polyvals.cols - 1) = graphMesh.y[i];
    }
    for (int j = 0; j < meshSizes.size(); j++) {
        auto mesh = meshGen(a, b, meshSizes[j], f);
        for (int i = 0; i <= graphSize; i++) {
            polyvals(i, j + 1) = newtonApprox(mesh, graphMesh.x[i]);
        }
    }
    ofstream polystream(dir + "polysmooth.bin", ios::binary);
    polyvals.write(polystream);
    polystream.close();

    matrix_t nvals(nhigh - nlow + 1);
    for (int i = nlow; i <= nhigh; i++) {
        auto mesh = meshGen(a, b, i, f);
        nvals[i - nlow] = maxDiv(mesh, f);
    }
    ofstream nstream(dir + "pointsmooth.bin", ios::binary);
    nvals.write(nstream);
    nstream.close();

    //функция с углом
    graphMesh = meshGen(a, b, graphSize, g);
    for (int i = 0; i <= graphSize; i++) {
        polyvals(i, 0) = graphMesh.x[i];
        polyvals(i, polyvals.cols - 1) = graphMesh.y[i];
    }

    for (int j = 0; j < meshSizes.size(); j++) {
        auto mesh = meshGen(a, b, meshSizes[j], g);
        for (int i = 0; i <= graphSize; i++) {
            polyvals(i, j + 1) = newtonApprox(mesh, graphMesh.x[i]);
        }
    }
    polystream.open(dir + "polycorner.bin", ios::binary);
    polyvals.write(polystream);
    polystream.close();

    for (int i = nlow; i <= nhigh; i++) {
        auto mesh = meshGen(a, b, i, g);
        nvals[i - nlow] = maxDiv(mesh, g);
        cout << "nval: " << nvals[i - nlow] << endl;
    }
    nstream.open(dir + "pointcorner.bin", ios::binary);
    nvals.write(nstream);
    nstream.close();
    return 0;
}

vector<double> divDiff(vector<pair<double, double>> vals) {
    int size = vals.size();
    matrix_t coef(size, size);
    for (int i = 0; i < size; i++) {
        coef(i, 0) = vals[i].second;
    }
    for (int i = 1; i < size; i++) {
        for (int j = 0; j < size - i; j++) {
            coef(j, i) = (coef(j, i - 1) - coef(j + 1, i - 1)) / (vals[j].first - vals[j + i].first);
        }
    }
    vector<double> ret(size);
    for (int i = 0; i < size; i++) {
        ret[i] = coef(0, i);
    }
    return ret;
}

double newtonApprox(vector<pair<double, double>> mesh, double x) {
    int size = mesh.size();
    vector<double> dd = divDiff(mesh);
    double res = dd[0];
    for (int i = 1; i < size; i++) {
        double prod = 1;
        for (int j = 0; j < i; j++) {
            prod *= x - mesh[j].first;
        }
        res += dd[i] * prod;
    }
    return res;
}

vector<pair<double, double>> meshGen(double a, double b, int n, double(f)(double)) {
    vector<pair<double, double>> mesh(n + 1);
    for (int i = 0; i <= n; i++) {
        mesh[i].first = a + (b - a) / n * i;
        mesh[i].second = f(mesh[i].first);
    }
    return mesh;
}

double maxDiv(vector<pair<double, double>> mesh, double(f)(double)) {
    int size = mesh.size();
    double d = 0;
    for (int i = 0; i < size - 1; i++) {
        double avg = (mesh[i].first + mesh[i + 1].first) / 2;
        double pval = newtonApprox(mesh, avg);
        d = max(d, abs(pval - f(avg)));
        cout << avg << ':' << pval << ':' << f(avg) << '\n';
    }
    return d;
}

mesh_t::mesh_t(int _size) : size(_size) {
    x.resize(_size);
    y.resize(_size);
}
