#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>


double Eps = 1e-5;
double eps = 1. / 2.;
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
        out << setprecision(10) << setw(25) << t << "|";
    }
    out << "\n" << "+";
    for (int i = 0; i < n; i++)
        out << "-------------------------+";
    out << "\n";
    while (flag) {
        for (int i = 0; i < n; i++) {
            in >> t;
            if (t == "/0") {
                flag = 0; break;
            }

            out << "|" << setprecision(10) << setw(25) << t;
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
    double n_r1cof_b_real = (D - A - B - C).real();
    double n_r1cof_b_imag = (D - A - B - C).imag();
    double n_r1cof_c_real = (D * D - A * D - B * D - C * D + A * B + A * C + B * C).real();
    double n_r1cof_c_imag = (D * D - A * D - B * D - C * D + A * B + A * C + B * C).imag();
    complex<double>E = second_root;
    double n_r2cof_b_real = (E + D - A - B - C).real();
    double n_r2cof_b_imag = (E + D - A - B - C).imag();
    double f;
    double x1 = x.real();
    double y1 = x.imag();
    if (!per1) {
        if (number_roots == 0)
            f = 2 * (-120 * pow(x1, 4) + 3 * pow(x1, 5) + 2 * pow(x1, 3) * (1165 - 50 * y1 + 3 * pow(y1, 2)) - 18 * pow(x1, 2) * (1443 - 138 * y1 + 8 * pow(y1, 2)) + x1 * (165240 - 27252 * y1 + 2402 * pow(y1, 2) - 100 * pow(y1, 3) + 3 * pow(y1, 4)) - 6 * (78300 - 19140 * y1 + 2283 * pow(y1, 2) - 138 * pow(y1, 3) + 4 * pow(y1, 4)));
        else if (number_roots == 1)
            f = 2 * (pow(n_r1cof_b_real, 2) * x1 + pow(n_r1cof_b_imag, 2) * x1 + 3 * n_r1cof_b_real * pow(x1, 2) + 2 * pow(x1, 3) + n_r1cof_c_real * (n_r1cof_b_real + 2 * x1) + 2 * n_r1cof_b_imag * x1 * y1 + (n_r1cof_b_real + 2 * x1) * pow(y1, 2) + n_r1cof_c_imag * (n_r1cof_b_imag + 2 * y1));
        else
            f = 2*(n_r2cof_b_real + x1);
        return (f);
    }
    else {
        if (number_roots == 0) 
            f = 2 * (-502200 + pow(x1, 3) * (828 - 96 * y1) + 178200 * y1 - 27918 * pow(y1, 2) + 2474 * pow(y1, 3) - 125 * pow(y1, 4) + 3 * pow(y1, 5) + pow(x1, 4) * (-25 + 3 * y1) + 2 * pow(x1, 2) * (-6813 + 1201 * y1 - 75 * pow(y1, 2) + 3 * pow(y1, 3)) - 12 * x1 * (-9570 + 2283 * y1 - 207 * pow(y1, 2) + 8 * pow(y1, 3)));
        else if (number_roots == 1)
            f = 2 * (n_r1cof_b_imag * pow(x1, 2) + n_r1cof_c_imag * (n_r1cof_b_real + 2 * x1) + pow(n_r1cof_b_real, 2) * y1 + pow(n_r1cof_b_imag, 2) * y1 + 2 * n_r1cof_b_real * x1 * y1 + 2 * pow(x1, 2) * y1 + 3 * n_r1cof_b_imag * pow(y1, 2) + 2 * pow(y1, 3) - n_r1cof_c_real * (n_r1cof_b_imag + 2 * y1));
        else
            f =  2 *(n_r2cof_b_imag+y1);
        
        return (f);
    }
}

complex<double>Grad(complex<double> x, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    complex<double>gr = complex<double>(DFunc(x, bool(0), number_roots, first_root, second_root), DFunc(x, bool(1), number_roots, first_root, second_root));
    return gr;
}


/*double GR(complex<double> x, bool per1, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
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
*/

complex<double> PoKordin(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
    out << "PoKordin" << "\n" << " " << "i" << " " << "Step" << " " << "(x,y)" << " " << "F(x,y)" << " " << "(x,y)_prew" << " " << "F_prew(x,y)" << endl;
    double step = 0.9;
    complex<double> p_xn = start;
    double Prew = Func(start, number_roots, first_root, second_root);
    

    complex<double> xn = p_xn + complex<double>(step,0.);
    double Curr = Func(xn, number_roots, first_root, second_root);
    out << i << " " << step << " " << xn << " " << Curr << " " << p_xn << " " << Prew << endl;
    while (abs(Curr - Prew) >eps) {
        xn = p_xn + complex<double>(step, 0.);
        Curr = Func(xn, number_roots, first_root, second_root);
        if (Prew < Curr) {
            step = -step;

            p_xn = xn;
            Prew = Curr;
            xn = p_xn + complex<double>(step, 0.);
            Curr = Func(xn, number_roots, first_root, second_root);
            i++;
            out << i << " " << step << " " << xn << " " << Curr << " " << p_xn << " " << Prew << endl;

        }

        while (Prew > Curr) {
            p_xn = xn;
            Prew = Curr;
            xn = p_xn + complex<double>(step, 0.);
            Curr = Func(xn, number_roots, first_root, second_root);
            i++;
            out << i << " " << step << " " << xn << " " << Curr << " " << p_xn << " " << Prew << endl;
        }
        xn = p_xn + complex<double>(0., step);
        Curr = Func(xn, number_roots, first_root, second_root);
        if (Prew < Curr) {
            step = -step;
            p_xn = xn;
            Prew = Curr;
            xn = p_xn + complex<double>(0., step);
            Curr = Func(xn, number_roots, first_root, second_root);
            i++;
            out << i << " " << step << " " << xn << " " << Curr << " " << p_xn << " " << Prew << endl;
        }
        while (Prew > Curr) {
            p_xn = xn;
            Prew = Curr;
            xn = p_xn + complex<double>(0., step);
            Curr = Func(xn, number_roots, first_root, second_root);
            i++;
            out << i << " " << step << " " << xn << " " << Curr << " " << p_xn << " " << Prew << endl;
        }
        
        step /= 2;
    }
    
    out << " /0 ";
    out << "Root at the point: " << xn << " Func value: " << " " << Curr << " /0 ";
    out << "Func was called " << coun << " times /0 \n";
    coun = 0;
    out.close();
    DrawTable(6);
    return (p_xn);
    
}

/*
complex<double> PoKordin2(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
    out << "PoKordin"<< number_roots+1 << "\n" << " " << "i" << " " << "(x,y)" << " " << "F(x,y)" << endl;
    
    double Prew = Func(start, number_roots, first_root, second_root);
    out <<i<< " " << start<< " " << Prew << endl;
    complex<double> xn = complex<double>(GR(start, bool(0) ,number_roots, first_root, second_root), start.imag());
    xn = complex<double>(xn.real(), GR(start, bool(1), number_roots, first_root, second_root));
    double Curr = Func(xn, number_roots, first_root, second_root);
    while (abs(Curr - Prew) > eps) {
        Prew = Curr;
        i++;
        out << i<<" " << xn << " " << Prew << endl;
        
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
*/

complex<double> GradWithConst(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
   complex<double>x = complex<double>(A, A);
    complex<double>y = complex<double>(B, B);
    double L = abs(Grad(x, number_roots, first_root, second_root) - Grad(y, number_roots, first_root, second_root))/abs(x-y);
    double a = 0.0001*pow(10,number_roots);
    
    double delta = 0.01;
    //double a = min((1 - eps) / (L),0.1);
    
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

complex<double> GradWithKnownStep(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
    //complex<double>x = complex<double>(A, A);
    //complex<double>y = complex<double>(B, B);
    //double L = abs(Grad(x, number_roots, first_root, second_root) - Grad(y, number_roots, first_root, second_root)) / abs(x - y);
    //double a = 0.0001*pow(10,number_roots);

    double delta = 0.01;
    double a = 1 / sqrt(i+1.);
    out << "GradWithKnownStep(a=" << a << ")" << "\n" << " " << "i" << " " << "(x,y)" << " " << "Grad(x,y)" << " " << "Abs(Grad(x,y))" << endl;

    complex<double>xn = start;
    complex<double>grad = Grad(xn, number_roots, first_root, second_root);
    while (abs(grad) >= delta) {

        out << i << " " << xn << " " << grad << " " << double(abs(grad)) << " " << endl;
        
        xn -= a * grad/ double(abs(grad));
        grad = Grad(xn, number_roots, first_root, second_root);
        i++;
        a = 1 / sqrt(i+1.);
    }
    out << i << " " << xn << " " << grad << " " << double(abs(grad)) << endl;
    out << " /0 ";
    out << "Root at the point: " << xn << " grad:" << " " << grad << "Func " << Func(xn) << " /0 ";
    out << "Func was called " << coun << " times /0 \n";
    coun = 0;
    out.close();
    DrawTable(4);
    return xn;
}

complex<double> GradWithCrushingStep(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
    
    //complex<double>x = complex<double>(A, A);
    //complex<double>y = complex<double>(B, B);
    //double L = abs(Grad(x, number_roots, first_root, second_root) - Grad(y, number_roots, first_root, second_root)) / abs(x - y);
    double a = 1, grad_abs, f_prew, f_xn;

    double delta = 0.01;
   // double a = 1 / sqrt(i + 1.);
    out << "GradWithCrushingStep" << "\n" << " " << "i" << " " << "a" << " " << "(x,y)" << " " << "Grad(x,y)" << " " << "Abs(Grad(x,y))" << " " << "F(xn)" << endl;

    complex<double>xn = start, prew;
    complex<double>grad = Grad(xn, number_roots, first_root, second_root);
    f_xn = Func(xn, number_roots, first_root, second_root);
    while (abs(grad) >= delta) {
        grad_abs = double(abs(grad));
        out << i << " " << a << " " << xn << " " << grad << " " << grad_abs << " " << f_xn << endl;
        prew = xn;
        xn -= a * grad / double(abs(grad));
        grad = Grad(xn, number_roots, first_root, second_root);
        i++;
        f_prew = f_xn;
        f_xn = Func(xn, number_roots, first_root, second_root);
        if (((f_prew - f_xn) <= -a * 0.9 * pow(grad_abs, 2)));
        else   a = a * 0.9;
    }
    out << i << " " << a << " " << xn << " " << grad << " " << double(abs(grad)) << " " << f_xn << endl;
    out << " /0 ";
    out << "Root at the point: " << xn << " grad:" << " " << grad << "Func " << Func(xn) << " /0 ";
    out << "Func was called " << coun << " times /0 \n";
    coun = 0;
    out.close();
    DrawTable(6);
    return xn;
}

complex<double> MNGS(complex<double> start, int number_roots = 0, complex<double> first_root = complex<double>(0, 0), complex<double> second_root = complex<double>(0, 0)) {
    int i = 0;
    ofstream out("tmptable.txt", ios_base::app);
    complex<double>xn = start;
    complex<double>grad = Grad(xn, number_roots, first_root, second_root);
    double func_prew;
    double func;
    double F;
    double step = 0.1;
    out << "MNGS" << "\n" << " " << "i" << " " << "(x,y)" << " " << "Grad(x,y)" << " " << "alpha" << " " << "F_prew((x,y)-alpha*grad)" << " " << "F((x,y)-alpha*grad)" << " " << "F(x,y)" << endl;
    double alpha,beta;
    bool flag = 1;
    while (abs(grad) >= eps) {
        
        alpha = 0.;
        if (flag)
            func_prew = Func(xn, number_roots, first_root, second_root);
        F = func_prew;
        beta = alpha + step;
        func = Func(xn - beta * grad, number_roots, first_root, second_root);
        out << i++ << " " << xn << " " << grad << " " << beta << " " << func_prew << " " << func << " " << F << endl;
        flag = 0;
        while (func_prew > func) {
            
            func_prew = func;
            alpha = beta;
            beta = alpha + step;
            func = Func(xn - beta * grad, number_roots, first_root, second_root);
            //i++;
            flag = 1;
            out << i++ << " " << xn << " " << grad << " " << beta << " " << func_prew << " " << func << " " << F << endl;
        }
        xn = xn - alpha * grad;
        F = func_prew;
        if (flag) {
            grad = Grad(xn, number_roots, first_root, second_root);
            alpha = 0;
            step = 0.1;
        }
        step /= 2;
    }
    out << i++ << " " << xn << " " << grad << " " << beta << " " << func_prew << " " << func << " " << F << endl;
    out << " /0 ";
    out << "Root at the point: " << xn << " grad:" << " " << grad << "Func " << Func(xn) << " /0 ";
    out << "Func was called " << coun << " times /0 \n";
    coun = 0;
    out.close();
    DrawTable(7);
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
    complex<double> second_root = PoKordin(first_root, 1, first_root);
    PoKordin(second_root, 2, first_root, second_root);

    
    first_root = GradWithConst(complex<double>(12., 12.), 0);
    second_root = GradWithConst(first_root, 1, first_root);
   GradWithConst(second_root, 2, first_root, second_root);

   first_root = GradWithKnownStep(complex<double>(12., 12.), 0);
   second_root = GradWithKnownStep(first_root, 1, first_root);
   GradWithKnownStep(second_root, 2, first_root, second_root);


   first_root = GradWithCrushingStep(complex<double>(12., 12.), 0);
   second_root = GradWithCrushingStep(first_root, 1, first_root);
   GradWithCrushingStep(second_root, 2, first_root, second_root);

   first_root = MNGS(complex<double>(0., 0.), 0);
   second_root = MNGS(first_root, 1, first_root);
   MNGS(second_root, 2, first_root, second_root);

	return 0;
}
