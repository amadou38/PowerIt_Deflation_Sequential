#include "headers.hpp"



int main(int argc, char *argv[])
{
	cout.precision(2);
    Matrix X(1,1);
    Vector x(1);
	Problem p = {p.A = X, p.v = x, .lambda = 0, p.Lmbd = x};

    int n = 17;   // taille de la matrice
    // p.A.resize(1,1);
    Tests(p, 3, n);
    Matrix A = p.A;
    Vector v = p.v;

    int m = 17;
    double tol = 1e-15;
    Matrix B(1,1);
    Vector Lambda(1);

    sortVect(p.Lmbd);
    
    cout << "Data present in the file: " << endl << endl;
    cout << "Matrix p.A (global): \n" << A << endl << endl;
    cout << "Initial vector p.v (global): \n" << v << endl << endl;
    cout << "\nExact Eigenvalues : " << endl;
    for (int i = 0; i < p.Lmbd.rows(); i++)
        cout << p.Lmbd(i) << "  ";
    cout << endl;
    
    auto t1 = std::chrono::system_clock::now();
    // Deflation method (with Puissance It)
    Deflation(p, B, m, tol);
    auto t2 = std::chrono::system_clock::now();

    chrono::duration<double> T = t2-t1;

    cout << "Approximate Eigenvalues (for Deflation): " << endl;
    for (int i = 0; i < m; i++)
        cout << p.Lmbd(i) << "  ";
    cout << "\nAnd associated eigenvectors (global):\n" << B << endl << endl;
    
    cout << "\n\nPuissance It & Deflation runtime: " << T.count() << " sec\n\n";
    
    int jj = 1;
    int maxiter = 100;
    Matrix Vm(1,1);
    Vector LmbdApp(1);
    TestERAM(p, LmbdApp, Vm, m, jj, tol, maxiter);

	return 0;
}