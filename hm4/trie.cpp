#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>

using namespace std;

struct TrieNode {
    map<string, TrieNode*> children;
    double coeff;
    bool is_end;
    
    TrieNode() : coeff(0.0), is_end(false) {}
};

class Polynomial {
private:
    TrieNode* root;
    
    void clearTrie(TrieNode* node) {
        if (!node) 
            return;
        for (auto& child : node->children) {
            clearTrie(child.second);
        }
        delete node;
    }
    
    TrieNode* copyTrie(TrieNode* node) {
        if (!node) 
            return nullptr;
        TrieNode* new_node = new TrieNode();
        new_node->coeff = node->coeff;
        new_node->is_end = node->is_end;
        for (auto& child : node->children) {
            new_node->children[child.first] = copyTrie(child.second);
        }
        return new_node;
    }
    
    map<string, int> parseMonomial(const string& monomial) {
        map<string, int> degrees;
        size_t pos = 0;
        
        while (pos < monomial.length()) {
            if (monomial[pos] == '*') {
                pos++;
                continue;
            }
            
            size_t start = pos;
            while (pos < monomial.length() && isalpha(monomial[pos])) {
                pos++;
            }
            
            if (start < pos) {
                string var_name = monomial.substr(start, pos - start);
                int power = 1;
                
                if (pos < monomial.length() && monomial[pos] == '^') {
                    pos++;
                    size_t pow_start = pos;
                    while (pos < monomial.length() && isdigit(monomial[pos])) {
                        pos++;
                    }
                    power = stoi(monomial.substr(pow_start, pos - pow_start));
                }
                
                degrees[var_name] = power;
            }
        }
        
        return degrees;
    }
    
    void insertMonomial(const map<string, int>& degrees, double coeff) {
        if (fabs(coeff) < 1e-10) 
            return;
        
        TrieNode* current = root;
        
        vector<pair<string, int>> sorted_degrees(degrees.begin(), degrees.end());
        sort(sorted_degrees.begin(), sorted_degrees.end());
        
        for (const auto& deg : sorted_degrees) {
            string key = deg.first + ":" + to_string(deg.second);
            if (current->children.find(key) == current->children.end()) {
                current->children[key] = new TrieNode();
            }
            current = current->children[key];
        }
        
        current->is_end = true;
        current->coeff += coeff;
        
        if (fabs(current->coeff) < 1e-10) {
            current->is_end = false;
            current->coeff = 0.0;
        }
    }
    
    void collectMonomials(TrieNode* node, vector<pair<map<string, int>, double>>& monomials,
                          map<string, int> current_degrees) const {
        if (!node) 
            return;
        
        if (node->is_end && fabs(node->coeff) > 1e-10) {
            monomials.push_back({current_degrees, node->coeff});
        }
        
        for (const auto& child : node->children) {
            size_t colon_pos = child.first.find(':');
            string var_name = child.first.substr(0, colon_pos);
            int power = stoi(child.first.substr(colon_pos + 1));
            
            map<string, int> new_degrees = current_degrees;
            new_degrees[var_name] = power;
            collectMonomials(child.second, monomials, new_degrees);
        }
    }
    
    string formatNumber(double num) const {
        if (fabs(num) < 1e-10) 
            return "";
        
        double int_part;
        double frac_part = modf(num, &int_part);
        
        if (fabs(frac_part) < 1e-10) {
            return to_string((long long)int_part);
        } else {
            string str = to_string(num);
            str.erase(str.find_last_not_of('0') + 1, string::npos);
            if (str.back() == '.') 
                str.pop_back();
            return str;
        }
    }
    
public:
    Polynomial() {
        root = new TrieNode();
    }
    
    Polynomial(const Polynomial& other) {
        root = copyTrie(other.root);
    }
    
    Polynomial& operator=(const Polynomial& other) {
        if (this != &other) {
            clearTrie(root);
            root = copyTrie(other.root);
        }
        return *this;
    }
    
    ~Polynomial() {
        clearTrie(root);
    }
    
    static Polynomial fromString(const string& str) {
        Polynomial result;
        string s = str;
        s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
        
        if (s.empty()) 
            return result;
        
        size_t pos = 0;
        double sign = 1.0;
        bool first_term = true;
        
        while (pos < s.length()) {
            if (s[pos] == '+') {
                sign = 1.0;
                pos++;
                first_term = false;
            } else if (s[pos] == '-') {
                sign = -1.0;
                pos++;
                first_term = false;
            } else if (first_term) {
                sign = 1.0;
                first_term = false;
            }
            
            double coeff = 1.0;
            bool has_number = false;
            size_t start = pos;
            
            if (pos < s.length() && (isdigit(s[pos]) || s[pos] == '.')) {
                has_number = true;
                while (pos < s.length() && (isdigit(s[pos]) || s[pos] == '.')) {
                    pos++;
                }
                coeff = stod(s.substr(start, pos - start));
            }
            
            string monomial_str;
            while (pos < s.length() && s[pos] != '+' && s[pos] != '-') {
                monomial_str += s[pos];
                pos++;
            }
            
            coeff *= sign;
            
            if (!monomial_str.empty()) {
                map<string, int> degrees = result.parseMonomial(monomial_str);
                result.insertMonomial(degrees, coeff);
            } else if (has_number) {
                map<string, int> empty_degrees;
                result.insertMonomial(empty_degrees, coeff);
            }
        }
        
        return result;
    }
    
    Polynomial operator+(const Polynomial& other) const {
        Polynomial result;
        vector<pair<map<string, int>, double>> monomials;
        
        collectMonomials(root, monomials, map<string, int>());
        for (const auto& mon : monomials) {
            result.insertMonomial(mon.first, mon.second);
        }
        
        monomials.clear();
        collectMonomials(other.root, monomials, map<string, int>());
        for (const auto& mon : monomials) {
            result.insertMonomial(mon.first, mon.second);
        }
        
        return result;
    }
    
    Polynomial operator-(const Polynomial& other) const {
        Polynomial result;
        vector<pair<map<string, int>, double>> monomials;
        
        collectMonomials(root, monomials, map<string, int>());
        for (const auto& mon : monomials) {
            result.insertMonomial(mon.first, mon.second);
        }
        
        monomials.clear();
        collectMonomials(other.root, monomials, map<string, int>());
        for (const auto& mon : monomials) {
            result.insertMonomial(mon.first, -mon.second);
        }
        
        return result;
    }
    
    Polynomial operator*(const Polynomial& other) const {
        Polynomial result;
        vector<pair<map<string, int>, double>> monomials1, monomials2;
        
        collectMonomials(root, monomials1, map<string, int>());
        collectMonomials(other.root, monomials2, map<string, int>());
        
        for (const auto& m1 : monomials1) {
            for (const auto& m2 : monomials2) {
                map<string, int> new_degrees = m1.first;
                for (const auto& deg : m2.first) {
                    new_degrees[deg.first] += deg.second;
                }
                double new_coeff = m1.second * m2.second;
                result.insertMonomial(new_degrees, new_coeff);
            }
        }
        
        return result;
    }
    
    vector<map<string, int>> getSupport() const {
        vector<map<string, int>> support;
        vector<pair<map<string, int>, double>> monomials;
        collectMonomials(root, monomials, map<string, int>());
        
        for (const auto& mon : monomials) {
            if (fabs(mon.second) > 1e-10) {
                support.push_back(mon.first);
            }
        }
        return support;
    }
    
    bool operator==(const Polynomial& other) const {
        vector<pair<map<string, int>, double>> m1, m2;
        collectMonomials(root, m1, map<string, int>());
        collectMonomials(other.root, m2, map<string, int>());
        
        if (m1.size() != m2.size()) 
            return false;
        
        for (size_t i = 0; i < m1.size(); i++) {
            if (m1[i].first != m2[i].first) 
                return false;
            if (fabs(m1[i].second - m2[i].second) > 1e-10) 
                return false;
        }
        
        return true;
    }
    
    bool operator!=(const Polynomial& other) const {
        return !(*this == other);
    }
    
    double evaluate(const map<string, double>& point) const {
        vector<pair<map<string, int>, double>> monomials;
        collectMonomials(root, monomials, map<string, int>());
        
        double result = 0.0;
        for (const auto& mon : monomials) {
            double term_value = mon.second;
            for (const auto& deg : mon.first) {
                auto it = point.find(deg.first);
                if (it != point.end()) {
                    term_value *= pow(it->second, deg.second);
                } else if (deg.second > 0) {
                    term_value = 0;
                    break;
                }
            }
            result += term_value;
        }
        
        return result;
    }
    
    int getHomogeneityDegree() const {
        vector<pair<map<string, int>, double>> monomials;
        collectMonomials(root, monomials, map<string, int>());
        
        if (monomials.empty()) 
            return -1;
        
        int first_degree = 0;
        for (const auto& deg : monomials[0].first) {
            first_degree += deg.second;
        }
        
        for (const auto& mon : monomials) {
            int current_degree = 0;
            for (const auto& deg : mon.first) {
                current_degree += deg.second;
            }
            if (current_degree != first_degree) {
                return -1;
            }
        }
        
        return first_degree;
    }
    
    pair<Polynomial, Polynomial> splitHomogeneous(int degree) const {
        Polynomial homogeneous;
        Polynomial remainder;
        
        vector<pair<map<string, int>, double>> monomials;
        collectMonomials(root, monomials, map<string, int>());
        
        for (const auto& mon : monomials) {
            int total_degree = 0;
            for (const auto& deg : mon.first) {
                total_degree += deg.second;
            }
            
            if (total_degree == degree) {
                homogeneous.insertMonomial(mon.first, mon.second);
            } else {
                remainder.insertMonomial(mon.first, mon.second);
            }
        }
        
        return {homogeneous, remainder};
    }
    
    string toString() const {
        vector<pair<map<string, int>, double>> monomials;
        collectMonomials(root, monomials, map<string, int>());
        
        if (monomials.empty()) 
            return "0";
        
        string result;
        bool first = true;
        
        sort(monomials.begin(), monomials.end(), 
             [](const pair<map<string, int>, double>& a, 
                const pair<map<string, int>, double>& b) {
                 int deg_a = 0, deg_b = 0;
                 for (const auto& d : a.first) 
                    deg_a += d.second;
                 for (const auto& d : b.first) 
                    deg_b += d.second;
                 if (deg_a != deg_b) 
                    return deg_a > deg_b;
                 return a.first.size() > b.first.size();
             });
        
        for (const auto& mon : monomials) {
            if (fabs(mon.second) < 1e-10) 
                continue;
            
            string coeff_str = formatNumber(mon.second);
            if (coeff_str.empty()) 
                continue;
            
            string monomial_str;
            if (!mon.first.empty()) {
                bool first_var = true;
                for (const auto& deg : mon.first) {
                    if (!first_var) 
                        monomial_str += "*";
                    monomial_str += deg.first;
                    if (deg.second > 1) {
                        monomial_str += "^" + to_string(deg.second);
                    }
                    first_var = false;
                }
            }
            
            if (!first) {
                if (mon.second > 0) {
                    result += " + ";
                } else {
                    result += " - ";
                    coeff_str = formatNumber(-mon.second);
                }
            } else {
                if (mon.second < 0) {
                    result += "-";
                    coeff_str = formatNumber(-mon.second);
                }
                first = false;
            }
            
            if (monomial_str.empty()) {
                result += coeff_str;
            } else {
                if (coeff_str != "1" || mon.first.empty()) {
                    result += coeff_str;
                    if (!monomial_str.empty()) 
                        result += "*";
                }
                result += monomial_str;
            }
        }
        
        return result;
    }
};

void printMenu() {
    cout << "\n=== МЕНЮ ===" << endl;
    cout << "a - Арифметические операции (сложение, вычитание, умножение)" << endl;
    cout << "b - Нахождение supp(f)" << endl;
    cout << "c - Проверка равенства двух многочленов" << endl;
    cout << "d - Вычисление значения многочлена в точке" << endl;
    cout << "e - Определение степени однородности" << endl;
    cout << "f - Разделение на однородную часть заданной степени" << endl;
    cout << "q - Выход" << endl;
    cout << "Выберите команду: ";
}

int main() {
    char choice;
    
    cout << "Программа для работы с многочленами от любого количества переменных" << endl;
    cout << "Поддерживаемые операции: +, -, *, ^ (через запись var^2)" << endl;
    cout << "Примеры ввода: 2*x^2*y + 3*x*y^2 - 5,  a*b*c + 2*a^2,  x^3 - 2*x*y + y^3" << endl;
    
    do {
        printMenu();
        cin >> choice;
        cin.ignore();
        
        switch (choice) {
            case 'a': {
                cout << "Введите первый многочлен: ";
                string p1_str;
                getline(cin, p1_str);
                Polynomial p1 = Polynomial::fromString(p1_str);
                
                cout << "Введите второй многочлен: ";
                string p2_str;
                getline(cin, p2_str);
                Polynomial p2 = Polynomial::fromString(p2_str);
                
                cout << "\nСумма: " << (p1 + p2).toString() << endl;
                cout << "Разность: " << (p1 - p2).toString() << endl;
                cout << "Произведение: " << (p1 * p2).toString() << endl;
                break;
            }
            
            case 'b': {
                cout << "Введите многочлен: ";
                string p_str;
                getline(cin, p_str);
                Polynomial p = Polynomial::fromString(p_str);
                
                vector<map<string, int>> support = p.getSupport();
                cout << "Supp(f) = { ";
                for (size_t i = 0; i < support.size(); i++) {
                    if (support[i].empty()) {
                        cout << "1";
                    } else {
                        bool first_var = true;
                        for (const auto& deg : support[i]) {
                            if (!first_var) 
                                cout << " * ";
                            cout << deg.first;
                            if (deg.second > 1) 
                                cout << "^" << deg.second;
                            first_var = false;
                        }
                    }
                    if (i < support.size() - 1) 
                        cout << ", ";
                }
                cout << " }" << endl;
                break;
            }
            
            case 'c': {
                cout << "Введите первый многочлен: ";
                string p1_str;
                getline(cin, p1_str);
                Polynomial p1 = Polynomial::fromString(p1_str);
                
                cout << "Введите второй многочлен: ";
                string p2_str;
                getline(cin, p2_str);
                Polynomial p2 = Polynomial::fromString(p2_str);
                
                if (p1 == p2) {
                    cout << "Многочлены равны" << endl;
                } else {
                    cout << "Многочлены не равны" << endl;
                }
                break;
            }
            
            case 'd': {
                cout << "Введите многочлен: ";
                string p_str;
                getline(cin, p_str);
                Polynomial p = Polynomial::fromString(p_str);
                
                cout << "Введите значения переменных (в формате: x=2, y=3, z=1): ";
                string point_str;
                getline(cin, point_str);
                
                map<string, double> point;
                stringstream ss(point_str);
                string token;
                while (getline(ss, token, ',')) {
                    size_t eq_pos = token.find('=');
                    if (eq_pos != string::npos) {
                        string var = token.substr(0, eq_pos);
                        var.erase(remove_if(var.begin(), var.end(), ::isspace), var.end());
                        double val = stod(token.substr(eq_pos + 1));
                        point[var] = val;
                    }
                }
                
                double result = p.evaluate(point);
                cout << "Значение многочлена = " << result << endl;
                break;
            }
            
            case 'e': {
                cout << "Введите многочлен: ";
                string p_str;
                getline(cin, p_str);
                Polynomial p = Polynomial::fromString(p_str);
                
                int degree = p.getHomogeneityDegree();
                if (degree == -1) {
                    cout << "Многочлен НЕ является однородным" << endl;
                } else {
                    cout << "Многочлен однородный степени " << degree << endl;
                }
                break;
            }
            
            case 'f': {
                cout << "Введите многочлен: ";
                string p_str;
                getline(cin, p_str);
                Polynomial p = Polynomial::fromString(p_str);
                
                int degree;
                cout << "Введите степень однородности: ";
                cin >> degree;
                cin.ignore();
                
                pair<Polynomial, Polynomial> result = p.splitHomogeneous(degree);
                cout << "Однородная часть степени " << degree << ": " << result.first.toString() << endl;
                cout << "Остаток: " << result.second.toString() << endl;
                break;
            }
            
            case 'q':
                cout << "Выход из программы..." << endl;
                break;
                
            default:
                cout << "Неверная команда! Попробуйте снова." << endl;
        }
        
    } while (choice != 'q');
    
    return 0;
}