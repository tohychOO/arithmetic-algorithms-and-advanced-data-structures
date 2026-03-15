#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

double f(double x) {
    return pow(x - 1, 3) * sin(M_PI * x) * (cos(2 * M_PI * x) - 1);
}

double df(double x) {
    double h = 1e-8;
    return (f(x + h) - f(x - h)) / (2 * h);
}

void newton(double x0, double eps) {
    cout << "\nМетод Ньютона для кратного корня\n";
    int iter = 0;
    double x = x0;
    double x_prev;

    do {
        x_prev = x;
        x = x_prev - f(x_prev) / df(x_prev);
        cout << "iter " << setw(2) << ++iter << ": x = " << x << ", f(x) = " << f(x) << "\n";
    } while (fabs(x - x_prev) > eps);

    cout << "Корень: " << x << ", итераций: " << iter << "\n";
}

int main() {
    cout << fixed << setprecision(8);
    cout << "Метод Ньютона для поиска корней уравнения (x-1)^3 sin(pi x) (cos(2pi x) - 1) = 0\n";

    double eps = 1e-6;

    newton(0.5, eps);
    newton(1.5, eps);
    newton(2.5, eps);

    return 0;
}