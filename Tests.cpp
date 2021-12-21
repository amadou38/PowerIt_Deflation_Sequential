#include "headers.hpp"


void TestBLAS(Problem& p)
{

    auto t1 = chrono::system_clock::now();
    double a = BLAS1(p.v, p.v);
    auto t2 = chrono::system_clock::now();
    Vector s = BLAS2(p.A, p.v);
    auto t3 = chrono::system_clock::now();
    // Vector ss = BLAS3(p.A, p.v);
    auto t4 = chrono::system_clock::now();

    // cout << "Matrix:\n" << p.A << endl;
    // cout << "Vector:\n" << p.v << endl;
    // cout << "a = " << a << endl << endl;
    cout << "s:\n" << s << endl << endl;
    // cout << "ss:\n" << ss << endl << endl;

    chrono::duration<double> T1 = t2-t1, T2 = t3-t2, T3 = t4-t3;
	//  cout << "\n\nBLAS1 runtime: " << T1.count() << " sec\n\nBLAS2 runtime: " << T2.count() << " sec\n\nBLAS3 runtime: " << T3.count() << " sec\n\n";
}

void TestERAM(Problem& p, Vector& LmbdApp, Matrix& Vm, int m, int jj, double tol, int maxiter)
{
    auto t1 = chrono::system_clock::now();
    ERAM(p, LmbdApp, Vm, m, jj, maxiter, tol);
    auto t2 = chrono::system_clock::now();
    chrono::duration<double> T = t2-t1;
    sortVect(p.Lmbd);
    // cout << "Data present in the file: " << endl << endl;
    // cout << "Matrix p.A: \n" << p.A << endl << endl;
    // cout << "Initial vector p.v: \n" << p.v << endl << endl;
    
    cout << "Exact Eigenvalues: " << endl;
    for (int i = 0; i < p.Lmbd.rows(); i++)
        cout << p.Lmbd(i) << "  ";
    cout << endl;
    cout << "Approximate Eigenvalues: " << endl;
    for (int i = 0; i < LmbdApp.rows(); i++)
        cout << LmbdApp(i) << "  ";
    cout << endl;
    cout << "Approximate Eigenvectors: \n" << Vm << endl;
    
    cout << "\n\nERAM runtime: " << T.count() << " sec\n\n";
    
}