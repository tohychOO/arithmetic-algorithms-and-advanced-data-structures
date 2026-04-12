#include "polynomials.h"
#include <iostream>

using namespace std;

int main() {
    vector<string> varOrder = {"x", "y", "z"};
    Ordering ord = LEX;
    
    auto f = Polynomial::fromString("x^2*y - z");
    auto g = Polynomial::fromString("x*y - 1");
    
    auto lm_f = f.getLeadingMonomial(ord, varOrder);
    auto lm_g = g.getLeadingMonomial(ord, varOrder);
    
    // gamma = LCM(lm_f, lm_g)
    auto gamma = Polynomial::lcmMonomial(lm_f, lm_g);
    
    auto S_fg = f.S(g, ord, varOrder);
    auto mdeg_S = S_fg.getMultiDegree(ord, varOrder);
    
    cout << "f = x^2*y - z" << endl;
    cout << "g = x*y - 1" << endl;
    cout << "lm(f) = (" << lm_f["x"] << "," << lm_f["y"] << "," << lm_f["z"] << ")" << endl;
    cout << "lm(g) = (" << lm_g["x"] << "," << lm_g["y"] << "," << lm_g["z"] << ")" << endl;
    cout << "gamma = (" << gamma["x"] << "," << gamma["y"] << "," << gamma["z"] << ")" << endl;
    cout << "S(f, g) = " << S_fg.toStringWithOrdering(varOrder, ord) << endl;
    cout << "multideg(S) = (" << mdeg_S[0] << "," << mdeg_S[1] << "," << mdeg_S[2] << ")" << endl;
    
    bool less = false;
    for (size_t i = 0; i < varOrder.size(); i++) {
        int s = (i < mdeg_S.size()) ? mdeg_S[i] : 0;
        int g_val = (varOrder[i] == "x") ? gamma["x"] :
                    (varOrder[i] == "y") ? gamma["y"] : gamma["z"];
        if (s < g_val) {
            less = true;
            break;
        }
        if (s > g_val) break;
    }
    
    cout << "multideg(S) < gamma: " << (less ? "да" : "нет") << endl;
    
    return 0;
}