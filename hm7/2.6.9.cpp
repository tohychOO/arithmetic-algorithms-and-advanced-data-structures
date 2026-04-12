#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    Ordering ord1 = GRLEX;
    Ordering ord2 = INVLEX;
    Ordering ord3 = LEX;
    
    // a
    auto f1 = Polynomial::fromString("x^2 - y");
    auto g1 = Polynomial::fromString("x^3 - z");
    
    vector<Polynomial> basis1 = {f1, g1};
    bool isGB1 = isGroebnerBasis(basis1, ord1, varOrder);
    cout << "a) Базис {x^2 - y, x^3 - z} при grlex-упорядочении" << (isGB1 ? "" : " не") << " является базисом Грёбнера: " << endl;

    // b
    auto f2 = Polynomial::fromString("x^2 - y");
    auto g2 = Polynomial::fromString("x^3 - z");
    
    vector<Polynomial> basis2 = {f2, g2};
    bool isGB2 = isGroebnerBasis(basis2, ord2, varOrder);
    cout << "b) Базис {x^2 - y, x^3 - z} при invlex-упорядочении" << (isGB2 ? "" : " не") << " является базисом Грёбнера: " << endl;

    // c
    auto f3 = Polynomial::fromString("x*y^2 - x*z + y");
    auto g3 = Polynomial::fromString("x*y - z^2");
    auto h3 = Polynomial::fromString("x - y*z^4");
    
    vector<Polynomial> basis3 = {f3, g3, h3};
    bool isGB3 = isGroebnerBasis(basis3, ord1, varOrder);
    cout << "c) Базис {x*y^2 - x*y - z^2, x - y*z^4} при lex-упорядочении" << (isGB3 ? "" : " не") << " является базисом Грёбнера: " << endl;

    return 0;
}