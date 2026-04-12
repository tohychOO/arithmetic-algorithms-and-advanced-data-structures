#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    Ordering ord = LEX;
    
    // a
    cout << "a) f = 4*x^2*z - 7*y^2, g = x*y*z^2 + 3*x*z^4" << endl;
    auto f1 = Polynomial::fromString("4*x^2*z - 7*y^2");
    auto g1 = Polynomial::fromString("x*y*z^2 + 3*x*z^4");
    
    auto S1 = f1.S(g1, ord, varOrder);
    cout << "   S(f, g) = " << S1.toStringWithOrdering(varOrder, ord) << endl;
    cout << endl;
    
    // b
    cout << "b) f = x^4*y - z^2, g = 3*x*z^2 - y" << endl;
    auto f2 = Polynomial::fromString("x^4*y - z^2");
    auto g2 = Polynomial::fromString("3*x*z^2 - y");

    auto S2 = f2.S(g2, ord, varOrder);
    cout << "   S(f, g) = " << S2.toStringWithOrdering(varOrder, ord) << endl;
    cout << endl;

    // c
    cout << "c) f = x^7*y^2*z + 2*i*x*y*z, g = 2*x^7*y^2*z + 4" << endl;
    auto f3 = Polynomial::fromString("x^7*y^2*z + 2*i*x*y*z");
    auto g3 = Polynomial::fromString("2*x^7*y^2*z + 4");
    
    auto S3 = f3.S(g3, ord, varOrder);
    cout << "   S(f, g) = " << S3.toStringWithOrdering(varOrder, ord) << endl;
    cout << endl;
    
    // d
    cout << "d) f = x*y + z^3, g = z^2 - 3*z" << endl;
    auto f4 = Polynomial::fromString("x*y + z^3");
    auto g4 = Polynomial::fromString("z^2 - 3*z");

    auto S4 = f4.S(g4, ord, varOrder);
    cout << "   S(f, g) = " << S4.toStringWithOrdering(varOrder, ord) << endl;
    
    return 0;
}