#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

string intToPolynomial(int value, int n) {
    if (value == 0) 
        return "0";
    
    vector<string> terms;
    
    for (int i = n - 1; i >= 0; i--) {
        if (value & (1 << i)) {
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

int polynomialToInt(const string& poly, int n) {
    int result = 0;
    
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
            result |= 1;
        else if (term == "x") 
            result |= (1 << 1);
        else if (term.substr(0, 2) == "x^") {
            int degree = stoi(term.substr(2));
            if (degree >= n) {
                cerr << "Ошибка: степень " << degree << " превышает допустимую" << endl;
                return 0;
            }
            result |= (1 << degree);
        } else {
            cerr << "Ошибка: непонятный терм '" << term << "'" << endl;
            return 0;
        }
    }
    
    return result;
}

string toBinary(int value, int n) {
    string result = "";
    for (int i = n - 1; i >= 0; i--) {
        if (value & (1 << i))
            result += "1";
        else
            result += "0";
    }
    return result;
}

int main() {
    cout << "GF(2^n) - преобразование между формами" << endl;
    cout << endl;
    
    int n;
    cout << "Введите степень поля n (2 <= n <= 30): ";
    cin >> n;
    
    if (n > 30) {
        cout << "Слишком большое n, возьмём n = 30" << endl;
        n = 30;
    }
    
    cout << "\nВыберите действие:" << endl;
    cout << "1. Число -> многочлен" << endl;
    cout << "2. Многочлен -> число" << endl;
    
    int choice;
    cin >> choice;
    cin.ignore();
    
    if (choice == 1) {
        int value;
        cout << "Введите число (0-" << (1 << n) - 1 << "): ";
        cin >> value;
        
        cout << "\nРезультат:" << endl;
        cout << "Число: " << value << endl;
        cout << "Двоичный вид: " << toBinary(value, n) << endl;
        cout << "Многочлен: " << intToPolynomial(value, n) << endl;
    } else if (choice == 2) {
        string poly;
        cout << "Введите многочлен (например: x^2 + x + 1): ";
        getline(cin, poly);
        
        int value = polynomialToInt(poly, n);
        
        cout << "\nРезультат:" << endl;
        cout << "Многочлен: " << poly << endl;
        cout << "Двоичный вид: " << toBinary(value, n) << endl;
        cout << "Число: " << value << endl;
    }
    
    return 0;
}