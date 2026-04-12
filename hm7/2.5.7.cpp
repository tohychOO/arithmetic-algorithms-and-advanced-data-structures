#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    Ordering ord = GRLEX;
    
    auto f1 = Polynomial::fromString("x^4*y^2 - z^5");
    auto f2 = Polynomial::fromString("x^3*y^3 - 1");
    auto f3 = Polynomial::fromString("x^2*y^4 - 2z");
    
    vector<Polynomial> basis = {f1, f2, f3};
    bool isGB = isGroebnerBasis(basis, ord, varOrder);
    cout << "Базис {x^4*y^2 - z^5, x^3*y^3 - 1, x^2*y^4 - 2z}" << (isGB ? "" : " не") << " является базисом Грёбнера: " << endl;

    return 0;
}