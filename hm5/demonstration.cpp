#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <polynomials.h>

using namespace std;

void printMenu() {
    cout << "\nДоступные команды" << endl;
    cout << "a - Арифметические операции (сложение, вычитание, умножение)" << endl;
    cout << "b - Нахождение supp(f)" << endl;
    cout << "c - Проверка равенства двух многочленов" << endl;
    cout << "d - Вычисление значения многочлена в точке" << endl;
    cout << "e - Определение степени однородности" << endl;
    cout << "f - Разделение на однородную часть заданной степени" << endl;
    cout << "g - Определение lm(f), lt(f), lc(f), multideg(f) с заданным упорядочением" << endl;
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
            
            case 'g': {
                cout << "Введите многочлен: ";
                string p_str;
                getline(cin, p_str);
                Polynomial p = Polynomial::fromString(p_str);
                
                cout << "Введите порядок переменных (через пробел, например: x y z): ";
                string order_line;
                getline(cin, order_line);
                stringstream ss_order(order_line);
                vector<string> varOrder;
                string var;
                while (ss_order >> var) {
                    varOrder.push_back(var);
                }
                
                cout << "Выберите упорядочение:" << endl;
                cout << "1 - lex" << endl;
                cout << "2 - grlex" << endl;
                cout << "3 - grevlex" << endl;
                cout << "4 - invlex" << endl;
                cout << "5 - rinvlex" << endl;
                cout << "Введите номер: ";
                int ord_choice;
                cin >> ord_choice;
                cin.ignore();
                
                Ordering ord;
                string ord_name;
                switch (ord_choice) {
                    case 1: ord = LEX; ord_name = "lex"; break;
                    case 2: ord = GRLEX; ord_name = "grlex"; break;
                    case 3: ord = GREVLEX; ord_name = "grevlex"; break;
                    case 4: ord = INVLEX; ord_name = "invlex"; break;
                    case 5: ord = RINVLEX; ord_name = "rinvlex"; break;
                    default: ord = LEX; ord_name = "lex"; break;
                }
                
                vector<pair<map<string, int>, double>> monomials = p.getAllMonomials();
                if (monomials.empty()) {
                    cout << "Многочлен равен 0" << endl;
                    break;
                }
                
                // Компаратор для поиска максимума (возвращает true, если a < b)
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
                
                // Находим старший элемент
                auto max_it = max_element(monomials.begin(), monomials.end(), comparator);
                map<string, int> leading_monom = max_it->first;
                double leading_coeff = max_it->second;
                
                // Функция для красивого вывода монома
                auto monomToString = [&](const map<string, int>& monom) -> string {
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
                };
                
                // Форматирование коэффициента (убираем .0 если целое)
                auto coeffToString = [](double c) -> string {
                    if (c == (int)c) return to_string((int)c);
                    string s = to_string(c);
                    s.erase(s.find_last_not_of('0') + 1, string::npos);
                    if (s.back() == '.') s.pop_back();
                    return s;
                };
                
                // Вывод результатов с информацией о многочлене и упорядочении
                cout << "\n=== Результаты для многочлена " << p.toString() 
                    << " с " << ord_name << "-упорядочением ===" << endl;
                
                // lc(f) — старший коэффициент
                cout << "lc(f) = " << coeffToString(leading_coeff) << endl;
                
                // lm(f) — старший моном
                cout << "lm(f) = " << monomToString(leading_monom) << endl;
                
                // lt(f) = lc(f) * lm(f)
                cout << "lt(f) = ";
                if (leading_coeff == 1 && !leading_monom.empty()) {
                    cout << monomToString(leading_monom);
                } else if (leading_coeff == -1 && !leading_monom.empty()) {
                    cout << "-" << monomToString(leading_monom);
                } else if (leading_monom.empty()) {
                    cout << coeffToString(leading_coeff);
                } else {
                    cout << coeffToString(leading_coeff) << "*" << monomToString(leading_monom);
                }
                cout << endl;
                
                // multideg(f) — вектор степеней старшего монома
                cout << "multideg(f) = (";
                for (size_t i = 0; i < varOrder.size(); i++) {
                    int deg = leading_monom.count(varOrder[i]) ? leading_monom.at(varOrder[i]) : 0;
                    cout << deg;
                    if (i < varOrder.size() - 1) cout << ", ";
                }
                cout << ")" << endl;
                
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