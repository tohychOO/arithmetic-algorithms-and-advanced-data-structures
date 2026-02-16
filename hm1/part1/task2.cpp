#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using Polynomial = vector<double>;

struct Solution {
    bool exists;
    vector<double> coefficients;
};

Polynomial inputPolynomial(const string& name) {
    cout << "Введите коэффициенты многочлена " << name 
         << " (от младшей степени к старшей, через пробел): ";
    
    string line;
    getline(cin, line);
    stringstream ss(line);
    
    Polynomial p;
    double coeff;
    while (ss >> coeff) {
        p.push_back(coeff);
    }
    
    return p;
}

void printPolynomial(const Polynomial& p, const string& name) {
    cout << name << "(x) = ";
    bool first = true;
    for (int i = 0; i < p.size(); i++) {
        if (abs(p[i]) < 1e-10) continue;
        
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

vector<double> gaussSolve(vector<vector<double>> A, vector<double> b) {
    int m = A.size();
    int n = A[0].size();
    
    vector<vector<double>> aug(m, vector<double>(n + 1));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            aug[i][j] = A[i][j];
        aug[i][n] = b[i];
    }
    
    for (int col = 0; col < min(m, n); col++) {
        int maxRow = col;
        for (int row = col + 1; row < m; row++) {
            if (abs(aug[row][col]) > abs(aug[maxRow][col]))
                maxRow = row;
        }
        
        if (abs(aug[maxRow][col]) < 1e-10)
            continue;
        
        swap(aug[col], aug[maxRow]);
        
        double pivot = aug[col][col];
        for (int j = col; j <= n; j++)
            aug[col][j] /= pivot;
        
        for (int row = col + 1; row < m; row++) {
            double factor = aug[row][col];
            for (int j = col; j <= n; j++)
                aug[row][j] -= factor * aug[col][j];
        }
    }
    
    vector<double> x(n, 0.0);
    
    for (int i = m - 1; i >= 0; i--) {
        bool nonZeroRow = false;
        int firstNonZero = -1;

        for (int j = 0; j < n; j++) {
            if (abs(aug[i][j]) > 1e-10) {
                nonZeroRow = true;
                firstNonZero = j;
                break;
            }
        }
        
        if (nonZeroRow) {
            x[firstNonZero] = aug[i][n];
            for (int j = firstNonZero + 1; j < n; j++)
                x[firstNonZero] -= aug[i][j] * x[j];
        } else {
            if (abs(aug[i][n]) > 1e-10)
                return vector<double>();
        }
    }
    
    return x;
}

Solution checkLinearSpan(const Polynomial& f, const vector<Polynomial>& g_list) {
    int max_deg = f.size() - 1;
    for (const auto& g : g_list) {
        if (g.size() - 1 > max_deg) {
            max_deg = g.size() - 1;
        }
    }
    
    int m = max_deg + 1;
    int n = g_list.size();
    
    vector<vector<double>> A(m, vector<double>(n, 0.0));
    vector<double> b(m, 0.0);
    
    for (int j = 0; j < n; j++) {
        const auto& g = g_list[j];
        for (int i = 0; i < g.size(); i++) {
            A[i][j] = g[i];
        }
    }
    
    for (int i = 0; i < f.size(); i++) {
        b[i] = f[i];
    }
    
    vector<double> coeffs = gaussSolve(A, b);
    
    if (coeffs.empty())
        return {false, {}};
    
    return {true, coeffs};
}

vector<vector<double>> binomialCoeffs(int max_n) {
    vector<vector<double>> C(max_n + 1, vector<double>(max_n + 1, 0));
    for (int n = 0; n <= max_n; n++) {
        C[n][0] = C[n][n] = 1;
        for (int k = 1; k < n; k++) {
            C[n][k] = C[n-1][k-1] + C[n-1][k];
        }
    }
    return C;
}

Polynomial powerXMinusA(int k, double a, const vector<vector<double>>& C) {
    Polynomial result(k + 1, 0);
    
    for (int j = 0; j <= k; j++) {
        double coeff = C[k][j] * pow(-a, k - j);
        result[j] = coeff;
    }
    
    return result;
}

vector<double> expandInPowerXMinusA(const Polynomial& f, double a) {
    int n = f.size() - 1;
    auto C = binomialCoeffs(n);
    
    vector<Polynomial> basis;
    for (int k = 0; k <= n; k++)
        basis.push_back(powerXMinusA(k, a, C));
    
    Solution sol = checkLinearSpan(f, basis);

    return sol.coefficients;
}

int main() {
    cout << "Разложение многочлена по степеням (x - a)" << endl;
    cout << endl;
    
    cout << "Введите коэффициенты многочлена f(x) (от младшей степени к старшей): ";
    string line;
    getline(cin, line);
    stringstream ss(line);
    
    Polynomial f;
    double coeff;
    while (ss >> coeff)
        f.push_back(coeff);
    
    double a;
    cout << "Введите число a: ";
    cin >> a;
    
    cout << endl;
    printPolynomial(f, "f");
    cout << "a = " << a << endl;
    cout << endl;
    
    vector<double> c = expandInPowerXMinusA(f, a);
    
    if (a > 0)
        cout << "Разложение по степеням (x - " << a << "):" << endl;
    else if (a == 0)
        cout << "Разложение по степеням x:" << endl;
    else
        cout << "Разложение по степеням (x + " << -a << "):" << endl;
    cout << "f(x) = ";
    
    bool first = true;
    for (int k = 0; k < c.size(); k++) {
        if (abs(c[k]) < 1e-10) continue;
        
        if (!first && c[k] > 0) 
            cout << " + ";
        else if (!first && c[k] < 0) 
            cout << " - ";
        else if (c[k] < 0) 
            cout << "- ";
        
        if (k == 0 || abs(c[k]) != 1) {
            cout << abs(c[k]);
        }
        
        if (k > 0) {
            if (a > 0)
                cout << "(x - " << a << ")";
            else if (a == 0)
                cout << "x";
            else
                cout << "(x + " << -a << ")";
            if (k > 1) 
                cout << "^" << k;
        }
        
        first = false;
    }
    if (first) 
        cout << "0";
    cout << endl;
}