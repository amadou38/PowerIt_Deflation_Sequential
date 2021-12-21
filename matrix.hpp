#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "headers.hpp"

using std::vector;
using std::tuple;
using std::ostream;
using std::istream;
class Matrix {
private:
    int grows;
    int gcols;
    vector<vector<double> > m_matrix;
public:
    Matrix(int, int);
    Matrix(int, int, double);
    Matrix(const char *);
    Matrix(const Matrix &);
    ~Matrix();
    
    // Matrix Operations
    Matrix operator+(Matrix &);
    Matrix operator-(Matrix &);
    Matrix operator*(Matrix &);
    Matrix transpose();
    
    // Scalar Operations
    Matrix operator+(double);
    Matrix operator-(double);
    Matrix operator*(double);
    Matrix operator/(double);
    
    // Aesthetic Methods
    double& operator()(const int &, const int &);
    void print() const;
    int rows() const;
    int cols() const;
    void conservativeResize(int, int);
    void resize(int, int);
    void random(int, int);
    friend ostream& operator << (ostream &, const Matrix &);
    friend istream& operator >> (istream &, Matrix &);
    
    // Power Iteration
    tuple<Matrix, double, int> powerIter(int, double);
    
    // Deflation
    Matrix deflation(Matrix &, double&);
};

class Vector {
private:
    int grows;
    vector<double> m_vector;
public:
    Vector(int);
    Vector(int, double);
    Vector(const Vector &);
    ~Vector();
    
    // Vector Operations
    Vector operator+(Vector &);
    Vector operator-(Vector &);
    Vector operator*(Vector &);
    Vector transpose();
    
    // Scalar Operations
    Vector operator+(double);
    Vector operator-(double);
    Vector operator*(double);
    Vector operator/(double);
    
    // Aesthetic Methods
    double& operator()(const int &);
    void print() const;
    int rows() const;
    void conservativeResize(int);
    void resize(int);
    void random(int);
    friend ostream& operator << (ostream &, const Vector &);
    friend istream& operator >> (istream &, Vector &);
    
};

#endif
