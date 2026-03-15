#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double phi_a(double x, double a, double b, double c) {
    return c + a * sin(x) * sin(x) + b * cos(x) * cos(x);
}

double dphi_a(double x, double a, double b) {
    return (a - b) * sin(2 * x);
}

void check_a() {
    double a, b, c, x0;
    int steps;

    cout << "\nφ(x) = c + a sin^2 x + b cos^2 x\n";
    cout << "Введите a, b, c: ";
    cin >> a >> b >> c;
    cout << "Введите начальное приближение x0: ";
    cin >> x0;
    cout << "Введите число итераций: ";
    cin >> steps;

    double B = a - b;
    double max_d = fabs(B);
    cout << "\n|a - b| = " << max_d << " → " << (max_d < 1 ? "сходимость" : "возможна расходимость") << "\n";

    double x = x0;
    for (int i = 0; i <= steps; ++i) {
        cout << "x[" << i << "] = " << x << "\n";
        x = phi_a(x, a, b, c);
    }
}

double phi_b(double x, double a, double b, double c) {
    return c + a * exp(-b * x * x);
}

double dphi_b(double x, double a, double b) {
    return -2 * a * b * x * exp(-b * x * x);
}

void check_b() {
    double a, b, c, x0;
    int steps;

    cout << "\nφ(x) = c + a e^{-b x^2}\n";
    cout << "Введите a, b > 0, c: ";
    cin >> a >> b >> c;
    cout << "Введите начальное приближение x0: ";
    cin >> x0;
    cout << "Введите число итераций: ";
    cin >> steps;

    double max_d = fabs(a) * sqrt(b) * sqrt(2.0) * exp(-0.5);
    cout << "\nmax |φ'(x)| ≈ " << max_d << " → " << (max_d < 1 ? "сходимость" : "возможна расходимость") << "\n";

    double x = x0;
    for (int i = 0; i <= steps; ++i) {
        cout << "x[" << i << "] = " << x << "\n";
        x = phi_b(x, a, b, c);
    }
}

int main() {
    cout << fixed << setprecision(6);
    cout << "Исследование сходимости метода простой итерации\n";

    int choice;
    do {
        cout << "\nВыберите пункт:\n";
        cout << "1) φ(x) = c + a sin^2 x + b cos^2 x\n";
        cout << "2) φ(x) = c + a e^{-b x^2}\n";
        cout << "0) Выход\n";
        cout << "Ваш выбор: ";
        cin >> choice;

        switch (choice) {
            case 1: 
                check_a(); 
                break;
            case 2: 
                check_b(); 
                break;
            case 0: 
                cout << "Выход.\n"; 
                break;
            default: 
                cout << "Неверный выбор.\n";
        }
    } while (choice != 0);

    return 0;
}