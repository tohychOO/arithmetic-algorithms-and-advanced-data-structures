#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder1 = {"x", "y"};
    vector<string> varOrder2 = {"y", "x"};
    vector<string> varOrder4 = {"x", "y", "z"};
    Ordering ord = LEX;
    
    // i
    auto f1 = Polynomial::fromString("x + y");
    auto g1 = Polynomial::fromString("y^2 - 1");
    
    vector<Polynomial> basis1 = {f1, g1};
    bool isGB1 = isGroebnerBasis(basis1, ord, varOrder1);
    cout << "  i) Базис {x + y, y^2 - 1} при x > y" << (isGB1 ? "" : " НЕ") << " является базисом Грёбнера: " << endl;

    // ii
    auto f2 = Polynomial::fromString("y + x");
    auto g2 = Polynomial::fromString("y^2 - 1");
    
    vector<Polynomial> basis2 = {f2, g2};
    bool isGB2 = isGroebnerBasis(basis2, ord, varOrder2);
    cout << " ii) Базис {y + x, y^2 - 1} при y > x" << (isGB2 ? "" : " НЕ") << " является базисом Грёбнера: " << endl;

    // iii
    auto f3 = Polynomial::fromString("x^2 + y^2 - 1");
    auto g3 = Polynomial::fromString("x*y - 1");
    auto h3 = Polynomial::fromString("x + y^3 - y");
    
    vector<Polynomial> basis3 = {f3, g3, h3};
    bool isGB3 = isGroebnerBasis(basis3, ord, varOrder1);
    cout << "iii) Базис {x^2 + y^2 - 1, x*y - 1, x + y^3 - y} при x > y" << (isGB3 ? "" : " НЕ") << " является базисом Грёбнера: " << endl;

    // iv
    auto f4 = Polynomial::fromString("x*y*z - 1");
    auto g4 = Polynomial::fromString("x - y");
    auto h4 = Polynomial::fromString("y*z^2 - 1");
    
    vector<Polynomial> basis4 = {f4, g4, h4};
    bool isGB4 = isGroebnerBasis(basis4, ord, varOrder4);
    cout << " iv) Базис {x*y*z - 1, x - y, y*z^2 - 1} при x > y > z" << (isGB4 ? "" : " НЕ") << " является базисом Грёбнера: " << endl;

    return 0;
}