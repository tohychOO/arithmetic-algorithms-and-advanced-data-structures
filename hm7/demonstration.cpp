#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    Ordering ord1 = GRLEX;
    Ordering ord2 = LEX;
    
    cout << "f = x^2, g = y^2" << endl;
    auto f1 = Polynomial::fromString("x^2");
    auto g1 = Polynomial::fromString("y^2");
    
    auto S1 = f1.S(g1, ord1, varOrder);
    cout << "S(f, g) = " << S1.toStringWithOrdering(varOrder, ord1) << endl;
    
    auto mdeg_f1 = f1.getMultiDegree(ord1, varOrder);
    auto mdeg_g1 = g1.getMultiDegree(ord1, varOrder);
    auto mdeg_fg1 = mdeg_f1;
    for (size_t i = 0; i < mdeg_fg1.size(); i++) {
        mdeg_fg1[i] += mdeg_g1[i];
    }
    cout << "multideg(f, g) = (";
    for (size_t i = 0; i < mdeg_fg1.size(); i++) {
        cout << mdeg_fg1[i];
        if (i + 1 < mdeg_fg1.size()) cout << ", ";
    }
    cout << ")" << endl;

    auto mdeg_S1 = S1.getMultiDegree(ord1, varOrder);
    cout << "multideg S(f, g) = (";
    for (size_t i = 0; i < mdeg_S1.size(); i++) {
        cout << mdeg_S1[i];
        if (i + 1 < mdeg_S1.size()) cout << ", ";
    }
    cout << ")" << endl;
    
    vector<Polynomial> basis1 = {f1, g1};
    bool isGB1 = isGroebnerBasis(basis1, ord1, varOrder);
    cout << "Базис является базисом Грёбнера: " << (isGB1 ? "да" : "нет") << endl;
    cout << endl;
    
    cout << "f = x^3*y^2 - x^2*y^3 + x, g = 3*x^4*y + y^2" << endl;
    auto f2 = Polynomial::fromString("x^3*y^2 - x^2*y^3 + x");
    auto g2 = Polynomial::fromString("3*x^4*y + y^2");

    auto S2 = f2.S(g2, ord2, varOrder);
    cout << "S(f, g) = " << S2.toStringWithOrdering(varOrder, ord2) << endl;
    
    auto mdeg_f2 = f2.getMultiDegree(ord2, varOrder);
    auto mdeg_g2 = g2.getMultiDegree(ord2, varOrder);
    auto mdeg_fg2 = mdeg_f2;
    for (size_t i = 0; i < mdeg_fg2.size(); i++) {
        mdeg_fg2[i] += mdeg_g2[i];
    }
    cout << "multideg(f, g) = (";
    for (size_t i = 0; i < mdeg_fg2.size(); i++) {
        cout << mdeg_fg2[i];
        if (i + 1 < mdeg_fg2.size()) cout << ", ";
    }
    cout << ")" << endl;
    
    auto mdeg_S2 = S2.getMultiDegree(ord2, varOrder);
    cout << "multideg S(f, g) = (";
    for (size_t i = 0; i < mdeg_S2.size(); i++) {
        cout << mdeg_S2[i];
        if (i + 1 < mdeg_S2.size()) cout << ", ";
    }
    cout << ")" << endl;
    
    vector<Polynomial> basis2 = {f2, g2};
    bool isGB2 = isGroebnerBasis(basis2, ord2, varOrder);
    cout << "Базис является базисом Грёбнера: " << (isGB2 ? "да" : "нет") << endl;
    
    return 0;
}