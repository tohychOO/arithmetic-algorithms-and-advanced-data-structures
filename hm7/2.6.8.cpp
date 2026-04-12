#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    Ordering ord = LEX;
    
    auto f1 = Polynomial::fromString("y - x^2");
    auto f2 = Polynomial::fromString("z - x^3");
    
    vector<Polynomial> basis = {f1, f2};
    bool isGB = isGroebnerBasis(basis, ord, varOrder);
    cout << "Базис {y - x^2, z - x^3}" << (isGB ? "" : " не") << " является базисом Грёбнера: " << endl;

    return 0;
}