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

int main() {
    cout << "Проверка принадлежности многочлена f(x) линейной оболочке, порожденной многочленами g1(x), g2(x), ... , gk(x)" << endl;
    cout << endl;

    Polynomial f = inputPolynomial("f(x)");
    int k;
    cout << "Введите количество многочленов g (базисных): ";
    cin >> k;
    cin.ignore();
    
    vector<Polynomial> g_list;
    for (int i = 0; i < k; i++) {
        Polynomial g = inputPolynomial("g" + to_string(i + 1));
        g_list.push_back(g);
    }
    
    cout << endl;
    printPolynomial(f, "f");
    for (int i = 0; i < g_list.size(); i++) {
        printPolynomial(g_list[i], "g" + to_string(i + 1));
    }
    cout << endl;
    
    Solution sol = checkLinearSpan(f, g_list);
    
    if (sol.exists) {
        cout << "f принадлежит линейной оболочке, порожденной многочленами g" << endl;
        cout << "Представление f через g:" << endl;
        cout << "f(x) = ";
        for (int i = 0; i < sol.coefficients.size(); i++) {
            if (i > 0) cout << " + ";
            cout << sol.coefficients[i] << "*g" << i + 1 << "(x)";
        }
        cout << endl;
    } else {
        cout << "f не принадлежит линейной оболочке, порожденной многочленами g" << endl;
        cout << "f нельзя представить через g" << endl;
    }
}