#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
using Polynomial = vector<double>;

void printPolynomial(const Polynomial& p, const string& name) {
    cout << name << "(x) = ";
    bool first = true;
    for (int i = 0; i < p.size(); i++) {
        if (abs(p[i]) < 1e-10) 
            continue;
        
        if (!first && p[i] > 0) 
            cout << " + ";
        else if (!first && p[i] < 0)
            cout << " - ";
        else if (p[i] < 0)
            cout << "- ";
        
        if (i == 0 || abs(p[i]) != 1) {
            cout << abs(p[i]);
        }
        
        if (i > 0) {
            cout << "x";
            if (i > 1) 
                cout << "^" << i;
        }
        
        first = false;
    }

    if (first) 
        cout << "0";
    cout << endl;
}

Polynomial invertSeries(const Polynomial& f, int n) {
    double a0 = f[0];
    Polynomial g(n + 1, 0.0);
    g[0] = 1.0 / a0;

    for (int i = 1; i < n; ++i) {
        double sum = 0.0;
        for (int k = 1; k <= i; ++k) {
            if (k < f.size())
                sum += f[k] * g[i - k];
        }
        g[i] = -sum / a0;
    }

    return g;
}

int main() {
    cout << setprecision(3);

    int n;
    cout << "Введите n > 0 (количество первых членов, которые необходимо найти): ";
    cin >> n;
    int M = n - 1;
    Polynomial exp_series(M + 1);
    double fact = 1.0;
    for (int i = 0; i <= M; ++i) {
        if (i > 0) 
            fact *= i;
        exp_series[i] = 1.0 / fact;
    }

    cout << "\nОбращение ряда для e^x (до x^" << n - 1 << "):\n";
    printPolynomial(exp_series, "f");
    auto inv_exp = invertSeries(exp_series, n);
    printPolynomial(inv_exp, "f^{-1}");

    Polynomial poly = {-1, -1, 1};
    cout << "\nОбращение ряда для x^2 - x - 1:\n";
    printPolynomial(poly, "f");
    auto inv_poly = invertSeries(poly, n);
    printPolynomial(inv_poly, "f^{-1}");

    return 0;
}