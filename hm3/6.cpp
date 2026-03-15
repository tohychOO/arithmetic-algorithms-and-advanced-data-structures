#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

int n;
double a;

double f_a(double x) {
    return sin(x) - 2*x*x + 0.5;
}
double df_a(double x) {
    return cos(x) - 4*x;
}

double f_b(double x, int n, double a) {
    double f = 1;
    for (int i = 0; i < n; ++i)
        f *= x;
    return f - a;
}

double f_b_wrapped(double x) {
    return f_b(x, n, a);
}

double df_b(double x, int n) {
    double f = n;
    for (int i = 0; i < n - 1; ++i)
        f *= x;
    return f;
}

double df_b_wrapped(double x) {
    return df_b(x, n);
}

double f_c(double x) {
    return sqrt(1 - x*x) - exp(x) + 0.1;
}
double df_c(double x) {
    return -x / sqrt(1 - x*x) - exp(x);
}

double f_d(double x) {
    return pow(x, 6) - 5*pow(x, 3) - 2;
}
double df_d(double x) {
    return 6*pow(x, 5) - 15*x*x;
}

double f_e(double x) {
    return log2(x) - 1.0/(1 + x*x);
}
double df_e(double x) {
    return 1.0/(x * log(2)) + (2*x)/((1 + x*x)*(1 + x*x));
}

double f_f(double x) {
    return sin(x/2) - 1;
}
double df_f(double x) {
    return 0.5 * cos(x/2);
}

double f_g(double x) {
    return log(x) - 1;
}
double df_g(double x) {
    return 1.0/x;
}

void dichotomy(double (*f)(double), double a, double b, double eps) {
    cout << "\nДихотомия\n";
    int iter = 0;
    double c;
    while ((b - a) > eps) {
        c = (a + b) / 2;
        cout << "iter " << setw(2) << ++iter << ": x = " << c << "\n";
        if (f(a) * f(c) <= 0)
            b = c;
        else
            a = c;
    }
    cout << "Корень: " << c << ", итераций: " << iter << "\n";
}

void Newton(double (*f)(double), double (*df)(double), double x0, double eps) {
    cout << "\nНьютон\n";
    int iter = 0;
    double x = x0;
    double x_prev;
    do {
        x_prev = x;
        x = x_prev - f(x_prev) / df(x_prev);
        cout << "iter " << setw(2) << ++iter << ": x = " << x << "\n";
    } while (fabs(x - x_prev) > eps);
    cout << "Корень: " << x << ", итераций: " << iter << "\n";
}

int main() {
    cout << fixed << setprecision(8);
    cout << "Нахождение корней методами дихотомии и Ньютона";

    cout << "Введите натуральное n (точность 10^(-n)): ";
    cin >> n;
    double eps = pow(10, -n);

    cout << "\na) sin(x) - 2x^2 + 0.5 = 0";
    dichotomy(f_a, -0.5, 0, eps);
    Newton(f_a, df_a, -0.25, eps);
    dichotomy(f_a, 0.5, 1, eps);
    Newton(f_a, df_a, 0.75, eps);

    cout << "\nb) x^n = a\n";
    cout << "Введите a > 0: ";
    cin >> a;
    if (a > 1) {
        if (n % 2 != 0) {
            dichotomy(f_b_wrapped, 1.0 + eps, a, eps);
            Newton(f_b_wrapped, df_b_wrapped, (1.0 + eps + a) / 2, eps);
        } else {
            dichotomy(f_b_wrapped, -a, -1.0 - eps, eps);
            Newton(f_b_wrapped, df_b_wrapped, (-a - 1.0 - eps) / 2, eps);
            dichotomy(f_b_wrapped, 1.0 + eps, a, eps);
            Newton(f_b_wrapped, df_b_wrapped, (1.0 + eps + a) / 2, eps);
        }
    } else if (a < 1) {
        if (n % 2 != 0) {
            dichotomy(f_b_wrapped, eps, 1.0 - eps, eps);
            Newton(f_b_wrapped, df_b_wrapped, 0.5, eps);
        } else {
            dichotomy(f_b_wrapped, -1.0 + eps, -eps, eps);
            Newton(f_b_wrapped, df_b_wrapped, -0.5, eps);
            dichotomy(f_b_wrapped, eps, 1.0 - eps, eps);
            Newton(f_b_wrapped, df_b_wrapped, 0.5, eps);
        } 
    } else {
        if (n % 2 != 0)
            cout << "Корень: 1" << endl;
        else 
            cout << "Корни: -1, 1" << endl;
    }

    cout << "\nc) sqrt(1-x^2) - e^x + 0.1 = 0";
    dichotomy(f_c, -1, -0.5, eps);
    Newton(f_c, df_c, -0.9, eps);
    dichotomy(f_c, 0, 0.5, eps);
    Newton(f_c, df_c, 0.25, eps);

    cout << "\nd) x^6 - 5x^3 - 2 = 0";
    dichotomy(f_d, -1, -0.5, eps);
    Newton(f_d, df_d, -0.75, eps);
    dichotomy(f_d, 1.5, 2, eps);
    Newton(f_d, df_d, 1.75, eps);

    cout << "\ne) log2(x) = 1/(1+x^2)";
    dichotomy(f_e, 1, 1.5, eps);
    Newton(f_e, df_e, 1.25, eps);

    cout << "\nf) sin(x/2) = 1";
    dichotomy(f_f, 3, 3.2, eps);
    Newton(f_f, df_f, 3.14, eps);

    cout << "\ng) ln(x) = 1";
    dichotomy(f_g, 2.7, 2.8, eps);
    Newton(f_g, df_g, 2.75, eps);

    return 0;
}