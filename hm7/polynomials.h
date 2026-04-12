#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>

using namespace std;

enum Ordering {
    LEX,
    GRLEX,
    GREVLEX,
    INVLEX,
    RINVLEX
};

int totalDegree(const map<string, int>& monom) {
    int deg = 0;
    for (const auto& p : monom)
        deg += p.second;
    return deg;
}

bool compareLex(const map<string, int>& a, const map<string, int>& b,
                const vector<string>& varOrder) {
    for (const string& var : varOrder) {
        int da = a.count(var) ? a.at(var) : 0;
        int db = b.count(var) ? b.at(var) : 0;
        if (da != db) 
            return da < db;
    }
    return false;
}

bool compareGrlex(const map<string, int>& a, const map<string, int>& b,
                  const vector<string>& varOrder) {
    int degA = totalDegree(a);
    int degB = totalDegree(b);
    if (degA != degB) 
        return degA < degB;
    return compareLex(a, b, varOrder);
}

bool compareGrevlex(const map<string, int>& a, const map<string, int>& b,
                    const vector<string>& varOrder) {
    int degA = totalDegree(a);
    int degB = totalDegree(b);
    if (degA != degB) 
        return degA < degB;
    for (auto it = varOrder.rbegin(); it != varOrder.rend(); ++it) {
        int da = a.count(*it) ? a.at(*it) : 0;
        int db = b.count(*it) ? b.at(*it) : 0;
        if (da != db) 
            return da > db;
    }
    return false;
}

bool compareInvlex(const map<string, int>& a, const map<string, int>& b,
                   const vector<string>& varOrder) {
    for (auto it = varOrder.rbegin(); it != varOrder.rend(); ++it) {
        int da = a.count(*it) ? a.at(*it) : 0;
        int db = b.count(*it) ? b.at(*it) : 0;
        if (da != db) 
            return da < db;
    }
    return false;
}

bool compareRinvlex(const map<string, int>& a, const map<string, int>& b,
                    const vector<string>& varOrder) {
    for (auto it = varOrder.rbegin(); it != varOrder.rend(); ++it) {
        int da = a.count(*it) ? a.at(*it) : 0;
        int db = b.count(*it) ? b.at(*it) : 0;
        if (da != db) 
            return da > db;
    }
    return !compareRinvlex(a, b, varOrder);
}

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
    
    vector<pair<map<string, int>, double>> getAllMonomials() const {
        vector<pair<map<string, int>, double>> monomials;
        collectMonomials(root, monomials, map<string, int>());
        return monomials;
    }

    void insertMonomial(const map<string, int>& degrees, double coeff) {
        if (fabs(coeff) < 1e-10) 
            return;
        
        TrieNode* current = root;
        
        vector<pair<string, int>> sorted_degrees(degrees.begin(), degrees.end());
        
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

    map<string, int> getLeadingMonomial(Ordering ord, const vector<string>& varOrder) const {
        auto monomials = getAllMonomials();
        if (monomials.empty()) return {};
        
        auto comparator = [&](const pair<map<string, int>, double>& a,
                            const pair<map<string, int>, double>& b) {
            switch (ord) {
                case LEX: return compareLex(a.first, b.first, varOrder);
                case GRLEX: return compareGrlex(a.first, b.first, varOrder);
                case GREVLEX: return compareGrevlex(a.first, b.first, varOrder);
                case INVLEX: return compareInvlex(a.first, b.first, varOrder);
                case RINVLEX: return compareRinvlex(a.first, b.first, varOrder);
                default: return compareLex(a.first, b.first, varOrder);
            }
        };
        
        auto max_it = max_element(monomials.begin(), monomials.end(), comparator);
        return max_it->first;
    }

    double getLeadingCoefficient(Ordering ord, const vector<string>& varOrder) const {
        auto monomials = getAllMonomials();
        if (monomials.empty()) return 0.0;
        
        auto comparator = [&](const pair<map<string, int>, double>& a,
                            const pair<map<string, int>, double>& b) {
            switch (ord) {
                case LEX: return compareLex(a.first, b.first, varOrder);
                case GRLEX: return compareGrlex(a.first, b.first, varOrder);
                case GREVLEX: return compareGrevlex(a.first, b.first, varOrder);
                case INVLEX: return compareInvlex(a.first, b.first, varOrder);
                case RINVLEX: return compareRinvlex(a.first, b.first, varOrder);
                default: return compareLex(a.first, b.first, varOrder);
            }
        };
        
        auto max_it = max_element(monomials.begin(), monomials.end(), comparator);
        return max_it->second;
    }

    pair<map<string, int>, double> getLeadingTerm(Ordering ord, const vector<string>& varOrder) const {
        auto monomials = getAllMonomials();
        if (monomials.empty()) return {{}, 0.0};
        
        auto comparator = [&](const pair<map<string, int>, double>& a,
                            const pair<map<string, int>, double>& b) {
            switch (ord) {
                case LEX: return compareLex(a.first, b.first, varOrder);
                case GRLEX: return compareGrlex(a.first, b.first, varOrder);
                case GREVLEX: return compareGrevlex(a.first, b.first, varOrder);
                case INVLEX: return compareInvlex(a.first, b.first, varOrder);
                case RINVLEX: return compareRinvlex(a.first, b.first, varOrder);
                default: return compareLex(a.first, b.first, varOrder);
            }
        };
        
        auto max_it = max_element(monomials.begin(), monomials.end(), comparator);
        return *max_it;
    }

    vector<int> getMultiDegree(Ordering ord, const vector<string>& varOrder) const {
        auto leading_monom = getLeadingMonomial(ord, varOrder);
        vector<int> result(varOrder.size(), 0);
        for (size_t i = 0; i < varOrder.size(); i++) {
            if (leading_monom.count(varOrder[i])) {
                result[i] = leading_monom.at(varOrder[i]);
            }
        }
        return result;
    }

    string monomialToString(const map<string, int>& monom, const vector<string>& varOrder) const {
        if (monom.empty()) return "1";
        string result;
        bool first = true;
        for (const string& v : varOrder) {
            int deg = monom.count(v) ? monom.at(v) : 0;
            if (deg > 0) {
                if (!first) result += "*";
                result += v;
                if (deg > 1) result += "^" + to_string(deg);
                first = false;
            }
        }
        return result.empty() ? "1" : result;
    }

    static string coeffToString(double c) {
        if (c == (int)c) return to_string((int)c);
        string s = to_string(c);
        s.erase(s.find_last_not_of('0') + 1, string::npos);
        if (s.back() == '.') s.pop_back();
        return s;
    }

    string toStringWithOrdering(const vector<string>& varOrder, Ordering ord) const {
        auto monomials = getAllMonomials();
        if (monomials.empty()) return "0";
        
        auto comparator = [&](const pair<map<string, int>, double>& a,
                              const pair<map<string, int>, double>& b) {
            switch (ord) {
                case LEX: return compareLex(b.first, a.first, varOrder);
                case GRLEX: return compareGrlex(b.first, a.first, varOrder);
                case GREVLEX: return compareGrevlex(b.first, a.first, varOrder);
                case INVLEX: return compareInvlex(b.first, a.first, varOrder);
                case RINVLEX: return compareRinvlex(b.first, a.first, varOrder);
                default: return compareLex(b.first, a.first, varOrder);
            }
        };
        
        sort(monomials.begin(), monomials.end(), comparator);
        
        string result;
        bool first = true;
        
        for (const auto& mon : monomials) {
            if (fabs(mon.second) < 1e-10) continue;
            
            string mon_str = monomialToString(mon.first, varOrder);
            string coeff_str = coeffToString(abs(mon.second));
            
            if (!first) {
                if (mon.second > 0) result += " + ";
                else result += " - ";
            } else {
                if (mon.second < 0) result += "-";
                first = false;
            }
            
            if (coeff_str == "1" && !mon.first.empty()) {
                result += mon_str;
            } else if (coeff_str == "1" && mon.first.empty()) {
                result += "1";
            } else if (mon.first.empty()) {
                result += coeff_str;
            } else {
                result += coeff_str + "*" + mon_str;
            }
        }
        
        return result;
    }

    static vector<pair<map<string, int>, double>> parseOrderedMonomials(const string& poly_str) {
        vector<pair<map<string, int>, double>> result;
        string s = poly_str;
        s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
        
        size_t pos = 0;
        double sign = 1.0;
        bool first = true;
        
        while (pos < s.length()) {
            if (s[pos] == '+') {
                sign = 1.0;
                pos++;
                first = false;
            } else if (s[pos] == '-') {
                sign = -1.0;
                pos++;
                first = false;
            } else if (first) {
                sign = 1.0;
                first = false;
            }
            
            double coeff = 1.0;
            bool has_number = false;
            size_t start = pos;
            
            if (pos < s.length() && (isdigit(s[pos]) || s[pos] == '.')) {
                has_number = true;
                while (pos < s.length() && (isdigit(s[pos]) || s[pos] == '.')) pos++;
                coeff = stod(s.substr(start, pos - start));
            }
            
            string monomial_str;
            while (pos < s.length() && s[pos] != '+' && s[pos] != '-') {
                monomial_str += s[pos];
                pos++;
            }
            
            coeff *= sign;
            
            if (!monomial_str.empty()) {
                map<string, int> degrees;
                size_t mpos = 0;
                while (mpos < monomial_str.length()) {
                    if (monomial_str[mpos] == '*') {
                        mpos++;
                        continue;
                    }
                    size_t vstart = mpos;
                    while (mpos < monomial_str.length() && isalpha(monomial_str[mpos])) mpos++;
                    string var = monomial_str.substr(vstart, mpos - vstart);
                    int power = 1;
                    if (mpos < monomial_str.length() && monomial_str[mpos] == '^') {
                        mpos++;
                        size_t pstart = mpos;
                        while (mpos < monomial_str.length() && isdigit(monomial_str[mpos])) mpos++;
                        power = stoi(monomial_str.substr(pstart, mpos - pstart));
                    }
                    degrees[var] = power;
                }
                result.push_back({degrees, coeff});
            } else if (has_number) {
                result.push_back({{}, coeff});
            }
        }
        
        return result;
    }

    static bool matchesOrdering(const vector<pair<map<string, int>, double>>& monomials,
                                const vector<string>& varOrder,
                                Ordering ord) {
        for (size_t i = 0; i < monomials.size() - 1; i++) {
            const auto& a = monomials[i].first;
            const auto& b = monomials[i + 1].first;
            
            bool a_less_b;
            switch (ord) {
                case LEX: a_less_b = compareLex(a, b, varOrder); break;
                case GRLEX: a_less_b = compareGrlex(a, b, varOrder); break;
                case GREVLEX: a_less_b = compareGrevlex(a, b, varOrder); break;
                default: a_less_b = compareLex(a, b, varOrder);
            }
            
            if (a_less_b) return false;
        }
        return true;
    }

    static string detectOrdering(const vector<pair<map<string, int>, double>>& monomials,
                                 const vector<string>& varOrder) {
        if (matchesOrdering(monomials, varOrder, LEX)) return "lex";
        if (matchesOrdering(monomials, varOrder, GRLEX)) return "grlex";
        if (matchesOrdering(monomials, varOrder, GREVLEX)) return "grevlex";
        return "не определено";
    }

    static bool divides(const map<string, int>& a, const map<string, int>& b) {
        for (const auto& p : a) {
            if (b.count(p.first) == 0 || b.at(p.first) < p.second)
                return false;
        }
        return true;
    }

    static map<string, int> divideMonomial(const map<string, int>& a, const map<string, int>& b) {
        map<string, int> result;
        for (const auto& p : b) {
            int deg_a = a.count(p.first) ? a.at(p.first) : 0;
            int deg_b = p.second;
            if (deg_b - deg_a > 0) {
                result[p.first] = deg_b - deg_a;
            }
        }
        return result;
    }

    struct DivisionResult;
    
    DivisionResult divide(const vector<Polynomial>& divisors, Ordering ord, 
                        const vector<string>& varOrder) const;

    static map<string, int> lcmMonomial(const map<string, int>& a, const map<string, int>& b) {
        map<string, int> result;
        for (const auto& p : a) result[p.first] = p.second;
        for (const auto& p : b) {
            if (result.count(p.first)) 
                result[p.first] = max(result[p.first], p.second);
            else 
                result[p.first] = p.second;
        }
        return result;
    }

    Polynomial S(const Polynomial& other, Ordering ord, const vector<string>& varOrder) const {
        auto lm1 = this->getLeadingMonomial(ord, varOrder);
        auto lm2 = other.getLeadingMonomial(ord, varOrder);
        auto lcm = lcmMonomial(lm1, lm2);
        
        auto lt1 = this->getLeadingTerm(ord, varOrder);
        auto lt2 = other.getLeadingTerm(ord, varOrder);
        
        auto factor1 = divideMonomial(lt1.first, lcm);
        auto factor2 = divideMonomial(lt2.first, lcm);
        
        Polynomial f1, f2;
        f1.insertMonomial(factor1, 1.0 / lt1.second);
        f2.insertMonomial(factor2, 1.0 / lt2.second);
        
        return (f1 * (*this)) - (f2 * other);
    }
};

struct Polynomial::DivisionResult {
    vector<Polynomial> quotients;
    Polynomial remainder;
};

Polynomial::DivisionResult Polynomial::divide(const vector<Polynomial>& divisors, Ordering ord, 
                                              const vector<string>& varOrder) const {
    int s = divisors.size();
    vector<Polynomial> quotients(s);
    Polynomial remainder;
    Polynomial p = *this;
    
    auto getLT = [&](const Polynomial& poly) -> pair<map<string, int>, double> {
        return poly.getLeadingTerm(ord, varOrder);
    };
    
    while (!p.getAllMonomials().empty()) {
        auto p_lt = getLT(p);
        bool divided = false;
        
        for (int i = 0; i < s; i++) {
            if (divisors[i].getAllMonomials().empty()) continue;
            
            auto f_lt = getLT(divisors[i]);
            
            if (divides(f_lt.first, p_lt.first)) {
                auto q_monom = divideMonomial(f_lt.first, p_lt.first);
                double q_coeff = p_lt.second / f_lt.second;
                
                Polynomial q_i;
                q_i.insertMonomial(q_monom, q_coeff);
                
                quotients[i] = quotients[i] + q_i;
                p = p - (q_i * divisors[i]);
                divided = true;
                break;
            }
        }
        
        if (!divided) {
            Polynomial lt_poly;
            lt_poly.insertMonomial(p_lt.first, p_lt.second);
            remainder = remainder + lt_poly;
            p = p - lt_poly;
        }
    }
    
    DivisionResult result;
    result.quotients = quotients;
    result.remainder = remainder;
    return result;
}

bool isGroebnerBasis(const vector<Polynomial>& basis, Ordering ord, const vector<string>& varOrder) {
    for (size_t i = 0; i < basis.size(); i++) {
        for (size_t j = i + 1; j < basis.size(); j++) {
            auto S_ij = basis[i].S(basis[j], ord, varOrder);
            auto divResult = S_ij.divide(basis, ord, varOrder);
            if (!divResult.remainder.getAllMonomials().empty()) {
                return false;
            }
        }
    }
    return true;
}

static void printResults(const Polynomial& p, const vector<string>& varOrder, Ordering ord, const string& ordName) {
    cout << "\n" << ordName << "-упорядочение" << endl;
    cout << "f = " << p.toStringWithOrdering(varOrder, ord) << endl;
    
    auto leading_term = p.getLeadingTerm(ord, varOrder);
    auto leading_monom = leading_term.first;
    double leading_coeff = leading_term.second;
    
    cout << "lc(f) = " << Polynomial::coeffToString(leading_coeff) << endl;
    cout << "lm(f) = " << p.monomialToString(leading_monom, varOrder) << endl;
    
    cout << "lt(f) = ";
    if (leading_coeff == 1 && !leading_monom.empty()) {
        cout << p.monomialToString(leading_monom, varOrder);
    } else if (leading_coeff == -1 && !leading_monom.empty()) {
        cout << "-" << p.monomialToString(leading_monom, varOrder);
    } else if (leading_monom.empty()) {
        cout << Polynomial::coeffToString(leading_coeff);
    } else {
        cout << Polynomial::coeffToString(leading_coeff) << "*" 
             << p.monomialToString(leading_monom, varOrder);
    }
    cout << endl;
    
    auto multideg = p.getMultiDegree(ord, varOrder);
    cout << "multideg(f) = (";
    for (size_t i = 0; i < multideg.size(); i++) {
        cout << multideg[i];
        if (i < multideg.size() - 1) cout << ", ";
    }
    cout << ")" << endl;
}