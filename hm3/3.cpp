#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

double f_a(double x) {
    return x*x*x + 3*x*x - 1;
}

double f_b(double x) {
    return x*x*x*x - x*x*x;
}

double f_c(double x) {
    return x*x - 3*x + 2;
}

double khi(double x, double alpha, double (*f)(double)) {
    return x + alpha * f(x);
}

void solve(double (*f)(double), const string& name) {
    cout << "\n" << name << "\n";

    double alpha, x0, eps;

    cout << "Введите alpha (< 0): ";
    cin >> alpha;

    cout << "Введите x0: ";
    cin >> x0;

    cout << "Введите точность: ";
    cin >> eps;

    cout << endl;
    cout << fixed << setprecision(8);

    double x_prev = x0;
    double x_curr, diff;
    int i = 0;

    do {
        ++i;
        x_curr = khi(x_prev, alpha, f);
        diff = abs(x_curr - x_prev);

        cout << "Шаг " << setw(2) << i << ": x = " << x_curr
             << " | f(x) = " << f(x_curr)
             << " | |Δ| = " << abs(x_curr - x_prev) << "\n";

        x_prev = x_curr;
    } while (diff > eps);

    cout << "-------------------------------------------------------------";
    cout << "\nПриближенный корень: " << x_curr << " (итераций: " << i << ")\n";
    return;
}

int main() {
    cout << "Метод простой итерации для трёх уравнений\n";
    cout << "-----------------------------------------\n";

    int choice;
    cout << "Выберите уравнение:\n";
    cout << "1) x^3 + 3x^2 - 1 = 0\n";
    cout << "2) x^4 - x^3 = 0\n";
    cout << "3) x^2 - 3x + 2 = 0\n";
    cout << "Ваш выбор: ";
    cin >> choice;
    cin.ignore();

    switch (choice) {
        case 1: 
            solve(f_a, "x^3 + 3x^2 - 1 = 0"); 
            break;
        case 2: 
            solve(f_b, "x^4 - x^3 = 0"); 
            break;
        case 3: 
            solve(f_c, "x^2 - 3x + 2 = 0"); 
            break;
        default: 
            cout << "Неверный выбор!\n";
    }

    return 0;
}