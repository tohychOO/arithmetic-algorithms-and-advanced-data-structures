#include <algorithm>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

uint64_t polynomialToInt(const string& poly, int n) {
    uint64_t result = 0;
    
    string s = poly;
    s.erase(remove(s.begin(), s.end(), ' '), s.end());
    
    if (s == "0") 
        return 0;
    
    stringstream ss(s);
    string term;
    
    while (getline(ss, term, '+')) {
        if (term.empty()) 
            continue;
        
        if (term == "1")
            result |= 1ULL;
        else if (term == "x") 
            result |= (1ULL << 1);
        else if (term.substr(0, 2) == "x^") {
            int degree = stoi(term.substr(2));
            if (degree >= n) {
                cerr << "Ошибка: степень " << degree << " превышает допустимую" << endl;
                return 0;
            }
            result |= (1ULL << degree);
        } else {
            cerr << "Ошибка: непонятный терм '" << term << "'" << endl;
            return 0;
        }
    }
    
    return result;
}

uint64_t multiplyPolynomials(uint64_t a, uint64_t b) {
    uint64_t result = 0;
    
    while (b != 0) {
        if (b & 1)
            result ^= a;
        
        a <<= 1;
        b >>= 1;
    }
    
    return result;
}

string intToPolynomial(uint64_t value, int n) {
    if (value == 0) 
        return "0";
    
    vector<string> terms;
    
    for (int i = n - 1; i >= 0; i--) {
        if (value & (1ULL << i)) {
            if (i == 0) {
                terms.push_back("1");
            } else if (i == 1) {
                terms.push_back("x");
            } else {
                terms.push_back("x^" + to_string(i));
            }
        }
    }
    
    string result;
    for (int i = 0; i < terms.size(); i++) {
        if (i > 0) 
            result += " + ";
        result += terms[i];
    }
    
    return result;
}

int main() {
    cout << "Умножение двоичных многочленов" << endl;
    cout << endl;
    cout << "Вводите многочлены в формате: x^2 + x + 1" << endl;

    cout << "Введите первый многочлен: ";
    string poly1;
    getline(cin, poly1);

    cout << "Введите второй многочлен: ";
    string poly2;
    getline(cin, poly2);
    
    int n = 32;
    uint64_t a = polynomialToInt(poly1, n);
    uint64_t b = polynomialToInt(poly2, n);
    
    cout << "\nИсходные данные:" << endl;
    cout << "a = " << poly1 << endl;
    cout << "b = " << poly2 << endl;
    
    uint64_t result = multiplyPolynomials(a, b);
    
    cout << "\nРезультат умножения:" << endl;
    cout << intToPolynomial(result, n * 2 + 1) << endl;
    
    return 0;
}