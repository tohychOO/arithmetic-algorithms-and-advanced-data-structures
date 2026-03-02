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
    uint64_t temp_a = a;
    uint64_t temp_b = b;
    
    while (temp_b != 0) {
        if (temp_b & 1)
            result ^= temp_a;
        
        temp_a <<= 1;
        temp_b >>= 1;
    }
    
    return result;
}

// ИСПРАВЛЕННАЯ ФУНКЦИЯ
uint64_t multiplyGF2n(uint64_t a, uint64_t b, uint64_t modulus, int n) {
    // Шаг 1: обычное умножение
    uint64_t product = multiplyPolynomials(a, b);
    
    // Шаг 2: редукция
    for (int i = 2*n - 2; i >= n; i--) {
        if (product & (1ULL << i)) {
            // Сдвигаем modulus так, чтобы его старший бит совпал с i
            uint64_t shifted_mod = modulus << (i - n);
            product ^= shifted_mod;
        }
    }
    
    return product;
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
    for (size_t i = 0; i < terms.size(); i++) {
        if (i > 0) 
            result += " + ";
        result += terms[i];
    }
    
    return result;
}

int main() {
    cout << "Умножение в GF(2^n)" << endl;
    cout << endl;
    
    int n;
    cout << "Введите степень поля n: ";
    cin >> n;
    cin.ignore();

    cout << "Введите первый элемент (многочлен степени < " << n << "):" << endl;
    string poly1;
    getline(cin, poly1);
    uint64_t a = polynomialToInt(poly1, n);
    
    cout << "Введите второй элемент (многочлен степени < " << n << "):" << endl;
    string poly2;
    getline(cin, poly2);
    uint64_t b = polynomialToInt(poly2, n);
    
    int num_mod;
    cout << "\nВведите количество неприводимых многочленов: ";
    cin >> num_mod;
    cin.ignore();
    
    vector<uint64_t> mod_values;
    vector<string> mod_strings;
    
    for (int i = 0; i < num_mod; i++) {
        cout << "Введите " << i+1 << "-й неприводимый многочлен степени " << n << ":" << endl;
        string mod_poly;
        getline(cin, mod_poly);
        
        uint64_t mod_val = polynomialToInt(mod_poly, n + 1);
        mod_values.push_back(mod_val);
        mod_strings.push_back(mod_poly);
    }

    for (int i = 0; i < num_mod; i++) {
        uint64_t result = multiplyGF2n(a, b, mod_values[i], n);
        
        cout << "\nНеприводимый многочлен " << i+1 << ": " << mod_strings[i] << " ---" << endl;
        cout << "Результат:" << endl;
        cout << intToPolynomial(result, n) << endl;
    }
    
    return 0;
}