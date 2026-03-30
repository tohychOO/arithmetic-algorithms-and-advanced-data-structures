#include <iostream>
#include <vector>
#include <string>
#include "polynomials.h"

using namespace std;

void printResults(const Polynomial& p, const vector<string>& varOrder, Ordering ord, const string& ordName) {
    cout << "\n" << ordName << "-упорядочение ===" << endl;
    
    auto leading_term = p.getLeadingTerm(ord, varOrder);
    auto leading_monom = leading_term.first;
    double leading_coeff = leading_term.second;
    
    cout << "lm(f) = " << p.monomialToString(leading_monom, varOrder) << endl;
    cout << "lt(f) = ";
    if (leading_coeff == 1 && !leading_monom.empty()) {
        cout << p.monomialToString(leading_monom, varOrder);
    } else if (leading_coeff == -1 && !leading_monom.empty()) {
        cout << "-" << p.monomialToString(leading_monom, varOrder);
    } else if (leading_monom.empty()) {
        cout << Polynomial::coeffToString(leading_coeff);
    } else {
        cout << Polynomial::coeffToString(leading_coeff) << "*" 
             << p.monomialToString(leading_monom, varOrder);
    }
    cout << endl;
    
    auto multideg = p.getMultiDegree(ord, varOrder);
    cout << "multideg(f) = (";
    for (size_t i = 0; i < multideg.size(); i++) {
        cout << multideg[i];
        if (i < multideg.size() - 1) cout << ", ";
    }
    cout << ")" << endl;
}

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    
    cout << "a) f = 2x + 3y + z + x^2 - z^2 + x^3" << endl;
    Polynomial p_a = Polynomial::fromString("2*x + 3*y + z + x^2 - z^2 + x^3");
    cout << "Многочлен: " << p_a.toString() << endl;
    
    printResults(p_a, varOrder, LEX, "lex");
    printResults(p_a, varOrder, GRLEX, "grlex");
    printResults(p_a, varOrder, GREVLEX, "grevlex");
    
    cout << "\n\nb) f = 2*x^2*y^8 - 3*x^5*y*z^4 + x*y^2*z^3 - x*y^4" << endl;
    Polynomial p_b = Polynomial::fromString("2*x^2*y^8 - 3*x^5*y*z^4 + x*y^2*z^3 - x*y^4");
    cout << "Многочлен: " << p_b.toString() << endl;
    
    printResults(p_b, varOrder, LEX, "lex");
    printResults(p_b, varOrder, GRLEX, "grlex");
    printResults(p_b, varOrder, GREVLEX, "grevlex");
    
    return 0;
}