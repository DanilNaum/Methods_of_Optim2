#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>


double Eps = 1e-5;
double eps = 1 / 2;
using namespace std;
int coun = 0;
double A = 0.;
double B = 13.;

void DrawTable(int n) {
    ofstream out("ans.txt", ios_base::app);
    string t;
    int flag = 1;
    ifstream in("tmptable.txt");
    in >> t;
    out << t << endl << "|";

    for (int i = 0; i < n; i++) {
        in >> t;
        out << setprecision(10) << setw(20) << t << "|";
    }
    out << "\n" << "+";
    for (int i = 0; i < n; i++)
        out << "--------------------+";
    out << "\n";
    while (flag) {
        for (int i = 0; i < n; i++) {
            in >> t;
            if (t == "/0") {
                flag = 0; break;
            }

            out << "|" << setprecision(10) << setw(20) << t;
        }
        if (flag)out << "|\n";
        else out << "\n";
    }
    flag = 1;
    while (flag) {
        in >> t;
        if (t == "/0") {
            flag = 0; break;
        }
        out << setprecision(3) << t << " ";
    }
    out << "\n\n";
    flag = 1;
    while (flag) {
        in >> t;
        if (t == "/0") {
            flag = 0; break;
        }
        out << setprecision(4) << t << " ";
    }
    out << "\n\n";

    in.close();
    out.close();
    ofstream out2("tmptable.txt");
    out2.close();
}

double Func(complex<double> x, int number_roots = 0, complex<double> first_root= complex<double>(0,0), complex<double> second_root=complex<double>(0,0)) {
	coun++;
    complex<double>A = complex<double>(6.0, 12.0);
    complex<double>B = complex<double>(12.0, 6.0);
    complex<double>C = complex<double>(6.0, 7.0);
    complex<double>D = first_root;
    complex<double>E = second_root;
    double f;
    if (number_roots == 0)
        f = abs((x - A) * (x - B) * (x - C));
    else if (number_roots == 1)
        f = abs(x * x + (D - A - B - C) * x + (D * D - A * D - B * D - C * D + A * B + A * C + B * C));
    else 
        f = abs(x + E + D - A - B - C);
    return (f);
}
double DFunc(complex<double> x, bool per1, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    coun++;

    complex<double>A = complex<double>(6.0, 12.0);
    complex<double>B = complex<double>(12.0, 6.0);
    complex<double>C = complex<double>(6.0, 7.0);
    complex<double>D = first_root;
    double dx = (D - A - B - C).real();
    double dy = (D - A - B - C).imag();
    double ddx = (D * D - A * D - B * D - C * D + A * B + A * C + B * C).real();
    double ddy = (D * D - A * D - B * D - C * D + A * B + A * C + B * C).imag();
    complex<double>E = second_root;
    double Ex = (E + D - A - B - C).real();
    double Ey = (E + D - A - B - C).imag();
    double f;
    double x1 = x.real();
    double y1 = x.imag();
    //if (per1) x1= complex<double>(x.real(),0.);
    //else x1 = complex<double>(0., x.imag());
    if (!per1) {
        if (number_roots == 0)
            f = 2 * (-120 * pow(x1, 4) + 3 * pow(x1, 5) + 2 * pow(x1, 3) * (1165 - 50 * y1 + 3 * pow(y1, 2)) - 18 * pow(x1, 2) * (1443 - 138 * y1 + 8 * pow(y1, 2)) + x1 * (165240 - 27252 * y1 + 2402 * pow(y1, 2) - 100 * pow(y1, 3) + 3 * pow(y1, 4)) - 6 * (78300 - 19140 * y1 + 2283 * pow(y1, 2) - 138 * pow(y1, 3) + 4 * pow(y1, 4)));
        else if (number_roots == 1)
            f = 2 * (pow(dx, 2) * x1 + pow(dy, 2) * x1 + 3 * dx * pow(x1, 2) + 2 * pow(x1, 3) + ddx * (dx + 2 * x1) + 2 * dy * x1 * y1 + (dx + 2 * x1) * pow(y1, 2) + ddy * (dy + 2 * y1));
        else
            f = 2*(Ex + x1);
        return (f);
    }
    else {
        if (number_roots == 0) 
            f = 2 * (-502200 + pow(x1, 3) * (828 - 96 * y1) + 178200 * y1 - 27918 * pow(y1, 2) + 2474 * pow(y1, 3) - 125 * pow(y1, 4) + 3 * pow(y1, 5) + pow(x1, 4) * (-25 + 3 * y1) + 2 * pow(x1, 2) * (-6813 + 1201 * y1 - 75 * pow(y1, 2) + 3 * pow(y1, 3)) - 12 * x1 * (-9570 + 2283 * y1 - 207 * pow(y1, 2) + 8 * pow(y1, 3)));
        else if (number_roots == 1)
            f = 2 * (dy * pow(x1, 2) + ddy * (dx + 2 * x1) + pow(dx, 2) * y1 + pow(dy, 2) * y1 + 2 * dx * x1 * y1 + 2 * pow(x1, 2) * y1 + 3 * dy * pow(y1, 2) + 2 * pow(y1, 3) - ddx * (dy + 2 * y1));
        else
            f =  2 *(Ey+y1);
        
        return (f);
    }
}
complex<double>Grad(complex<double> x, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    complex<double>gr = complex<double>(DFunc(x, bool(0), number_roots, first_root, second_root), DFunc(x, bool(1), number_roots, first_root, second_root));
    return gr;
}


double GR(complex<double> x, bool per1, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    double a = A;
    double b = B;
    double sqrt5 = sqrt(5);
    double c = (3 - sqrt5) / 2 * (b - a) + a;
    double d = (sqrt5 - 1) / 2 * (b - a) + a;
    
    double fc;
    double fd;
    if (per1) {fc = Func(complex<double>(x.real(), c), number_roots, first_root, second_root); fd = Func(complex<double>(x.real(), d), number_roots, first_root, second_root);}
    else {fc = Func(complex<double>(c, x.imag()), number_roots, first_root, second_root); fd = Func(complex<double>(d, x.imag()), number_roots, first_root, second_root);}
    
    while ((b - a) / 2 > Eps){
        if (fc <= fd) {
            b = d;
            d = c;
            fd = fc;
            c = (3 - sqrt(5)) * (b - a) / 2 + a;
            if (per1) { fc = Func(complex<double>(x.real(), c), number_roots, first_root, second_root); }
            else { fc = Func(complex<double>(c, x.imag()), number_roots, first_root, second_root);}


        }
        else {
            a = c;
            c = d;
            fc = fd;
            d = (sqrt(5) - 1) * (b - a) / 2 + a;
            if (per1) { fd = Func(complex<double>(x.real(), d), number_roots, first_root, second_root); }
            else { fd = Func(complex<double>(d, x.imag()), number_roots, first_root, second_root); }

        }
    }
    return((a+b)/2);
}

complex<double> PoKordin(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
    out << "PoKordin"<< number_roots+1 << "\n" << " " << "i" << " " << "(x,y)" << " " << "F(x,y)" << endl;
    //out << "Passive_Search" << "\n" << " " << "i" << " " << "x" << " " << "y" << endl;
    double Prew = Func(start, number_roots, first_root, second_root);
    out <<i<< " " << start<< " " << Prew << endl;
    complex<double> xn = complex<double>(GR(start, bool(0) ,number_roots, first_root, second_root), start.imag());
    xn = complex<double>(xn.real(), GR(start, bool(1), number_roots, first_root, second_root));
    double Curr = Func(xn, number_roots, first_root, second_root);
    while (abs(Curr - Prew) > eps) {
        Prew = Curr;
        i++;
        out << i<<" " << xn << " " << Prew << endl;
        //xn = complex<double>(GR(xn, bool(0), number_roots, first_root, second_root), GR(xn, bool(1), number_roots, first_root, second_root));
        xn = complex<double>(GR(xn, bool(0), number_roots, first_root, second_root), xn.imag());
        xn = complex<double>(xn.real(), GR(xn, bool(1), number_roots, first_root, second_root));
        Curr = Func(xn, number_roots, first_root, second_root);
    }
    out << " /0 ";
    out << "Root at the point: " << xn << " Func value: " << " " << Curr << " /0 ";
    out << "Func was called " << coun << " times /0 \n";
    coun = 0;
    out.close();
    DrawTable(3);
    return xn;
}

complex<double> GradWithConst(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
   complex<double>x = complex<double>(A, A);
    complex<double>y = complex<double>(B, B);
    double L = abs(Grad(x, number_roots, first_root, second_root)- Grad(y, number_roots, first_root, second_root))/abs(x-y);
    //double a = 0.0001*pow(10,number_roots);
    
    double delta = 0.01;
    double a = min((1 - eps) / (L),0.1);
    out << "GradWithConst(a=" << a << ")" << "\n" << " " << "i" << " " << "(x,y)" << " " << "Grad(x,y)" << " " << "Abs(Grad(x,y))" << endl;

    complex<double>xn = start;
    complex<double>grad = Grad(xn, number_roots, first_root, second_root);
    while(abs(grad)>= delta){
        
        out <<i<<" " << xn<<" " << grad <<" " << double(abs(grad)) << endl;
        xn -= a * grad;
        grad = Grad(xn, number_roots, first_root, second_root);
        i++;
    }
    out << i << " " << xn << " " << grad << " " << double(abs(grad)) << endl;
    out << " /0 ";
    out << "Root at the point: "  << xn << " grad:"  << " " << grad<< "Func "<< Func(xn) << " /0 ";
    out << "Func was called " << coun << " times /0 \n";
    coun = 0;
    out.close();
    DrawTable(4);
    return xn;
}
void clear_out(string file) {
    ofstream out(file);
    out.close();
    
    return;
}
int main() {
    clear_out("ans.txt");
    clear_out("tmptable.txt");
    complex<double> first_root = PoKordin(complex<double>(0.,0.),0);
    complex<double> second_root = PoKordin(complex<double>(0., 0.), 1, first_root);
    PoKordin(second_root, 2, first_root, second_root);

    
    first_root = GradWithConst(complex<double>(12., 12.), 0);
    second_root = GradWithConst(first_root, 1, first_root);
   GradWithConst(second_root, 2, first_root, second_root);


	return 0;
}
