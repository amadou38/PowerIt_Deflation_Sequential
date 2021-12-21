#include "headers.hpp"

//using std::ostream;
// Constructor for Any Matrix
Matrix::Matrix(int rowSize, int colSize){
    grows = rowSize;
    gcols = colSize;
    m_matrix.resize(rowSize);
    double init = 0.0;
    for (int i = 0; i < m_matrix.size(); i++)
    {
        m_matrix[i].resize(colSize, init);
    }
}

Matrix::Matrix(int rowSize, int colSize, double init){
    grows = rowSize;
    gcols = colSize;
    m_matrix.resize(rowSize);
    for (int i = 0; i < m_matrix.size(); i++)
    {
        m_matrix[i].resize(colSize, init);
    }
}

void Matrix::conservativeResize(int rowSize, int colSize)
{
    Matrix B(grows,gcols);
    for (int i = 0; i < grows; ++i)
        for (int j = 0; j < gcols; ++j)
            B(i,j) = m_matrix[i][j];
    grows = rowSize;
    gcols = colSize;
    m_matrix.resize(rowSize);
    double init = 0.0;
    for (int i = 0; i < m_matrix.size(); i++)
        m_matrix[i].resize(colSize, init);

    for (int i = 0; i < grows; ++i)
        for (int j = 0; j < gcols; ++j)
            this->m_matrix[i][j] = B(i,j);
}

void Matrix::resize(int rowSize, int colSize)
{
    grows = rowSize;
    gcols = colSize;
    m_matrix.resize(rowSize);
    double init = 0.0;
    for (int i = 0; i < m_matrix.size(); i++)
        m_matrix[i].resize(colSize, init);
}

void Matrix::random(int rowSize, int colSize)
{
    grows = rowSize;
    gcols = colSize;
    Matric B = Matric::Random(grows,gcols);
    m_matrix.resize(rowSize);
    double init = 0.0;
    for (int i = 0; i < m_matrix.size(); i++)
        m_matrix[i].resize(colSize, init);

    for (int i = 0; i < grows; ++i)
        for (int j = 0; j < gcols; ++j)
            this->m_matrix[i][j] = B(i,j);
}

// Constructor for Given Matrix
Matrix::Matrix(const char * fileName){
    ifstream file_A(fileName); // input file stream to open the file A.txt

    // Task 1
    // Keeps track of the Column and row sizes
    int colSize = 0;
    int rowSize = 0;
    
    // read it as a vector
    string line_A;
    int idx = 0;
    double element_A;
    double *vector_A = nullptr;
    
    if (file_A.is_open() && file_A.good())
    {
        // cout << "File A.txt is open. \n";
        while (getline(file_A, line_A))
        {
            rowSize += 1;
            stringstream stream_A(line_A);
            colSize = 0;
            while (1)
            {
                stream_A >> element_A;
                if (!stream_A)
                    break;
                colSize += 1;
                double *tempArr = new double[idx + 1];
                copy(vector_A, vector_A + idx, tempArr);
                tempArr[idx] = element_A;
                vector_A = tempArr;
                idx += 1;
            }
        }
    }
    else
    {
        cout << " Failed to open. \n";
    }
    
    int j;
    idx = 0;
    m_matrix.resize(rowSize);
    for (int i = 0; i < m_matrix.size(); i++) {
        m_matrix[i].resize(colSize);
    }
    for (int i = 0; i < rowSize; i++)
    {
        for (j = 0; j < colSize; j++)
        {
            this->m_matrix[i][j] = vector_A[idx];
            idx++;
        }
    }
    gcols = colSize;
    grows = rowSize;
    delete [] vector_A; // Tying up loose ends
    

}

// Copy Constructor
Matrix::Matrix(const Matrix &B)
{
    this->gcols = B.cols();
    this->grows = B.rows();
    this->m_matrix = B.m_matrix;
    
}

Matrix::~Matrix(){

}

// Addition of Two Matrices
Matrix Matrix::operator+(Matrix &B){
    Matrix sum(grows,B.cols(), 0.0);
    int i,j;
    for (i = 0; i < grows; i++)
    {
        for (j = 0; j < gcols; j++)
        {
            sum(i,j) = this->m_matrix[i][j] + B(i,j);
        }
    }
    return sum;
}

// Subtraction of Two Matrices
Matrix Matrix::operator-(Matrix & B){
    Matrix diff(grows,B.cols(), 0.0);
    int i,j;
    for (i = 0; i < grows; i++)
    {
        for (j = 0; j < gcols; j++)
        {
            diff(i,j) = this->m_matrix[i][j] - B(i,j);
        }
    }
    
    return diff;
}

// Multiplication of Two Matrices
Matrix Matrix::operator*(Matrix & B){
    Matrix multip(grows,B.cols(),0.0);
    if(gcols == B.rows())
    {
        int i,j,k;
        double temp = 0.0;
        for (i = 0; i < grows; i++)
        {
            for (j = 0; j < B.cols(); j++)
            {
                temp = 0.0;
                for (k = 0; k < gcols; k++)
                {
                    temp += m_matrix[i][k] * B(k,j);
                }
                multip(i,j) = temp;
                //cout << multip(i,j) << " ";
            }
            //cout << endl;
        }
        return multip;
    }
    else
    {
        return "Error";
    }
}

// Scalar Addition
Matrix Matrix::operator+(double scalar){
    Matrix result(grows,gcols,0.0);
    int i,j;
    for (i = 0; i < grows; i++)
    {
        for (j = 0; j < gcols; j++)
        {
            result(i,j) = this->m_matrix[i][j] + scalar;
        }
    }
    return result;
}

// Scalar Subraction
Matrix Matrix::operator-(double scalar){
    Matrix result(grows,gcols,0.0);
    int i,j;
    for (i = 0; i < grows; i++)
    {
        for (j = 0; j < gcols; j++)
        {
            result(i,j) = this->m_matrix[i][j] - scalar;
        }
    }
    return result;
}

// Scalar Multiplication
Matrix Matrix::operator*(double scalar){
    Matrix result(grows,gcols,0.0);
    int i,j;
    for (i = 0; i < grows; i++)
    {
        for (j = 0; j < gcols; j++)
        {
            result(i,j) = this->m_matrix[i][j] * scalar;
        }
    }
    return result;
}

// Scalar Division
Matrix Matrix::operator/(double scalar){
    Matrix result(grows,gcols,0.0);
    int i,j;
    for (i = 0; i < grows; i++)
    {
        for (j = 0; j < gcols; j++)
        {
            result(i,j) = this->m_matrix[i][j] / scalar;
        }
    }
    return result;
}


// Returns value of given location when asked in the form A(x,y)
double& Matrix::operator()(const int &rowNo, const int & colNo)
{
    return this->m_matrix[rowNo][colNo];
}

// No brainer - returns row #
int Matrix::rows() const
{
    return this->grows;
}

// returns col #
int Matrix::cols() const
{
    return this->gcols;
}

// Take any given matrices transpose and returns another matrix
Matrix Matrix::transpose()
{
    Matrix Transpose(gcols,grows,0.0);
    for (int i = 0; i < gcols; i++)
    {
        for (int j = 0; j < grows; j++) {
            Transpose(i,j) = this->m_matrix[j][i];
        }
    }
    return Transpose;
}

// Prints the matrix beautifully
size_t number_of_digits(double n) {
    std::ostringstream strs;

    strs << n;
    return strs.str().size();
}
ostream& operator << (ostream &output, const Matrix &A)
{

    int grows = A.grows, gcols = A.gcols;
    size_t max_len_per_column[grows];

    for (size_t j = 0; j < gcols; ++j) {
        size_t max_len {};

        for (size_t i = 0; i < grows; ++i)
            if (const auto num_length {number_of_digits(A.m_matrix[i][j])}; num_length > max_len)
                max_len = num_length;

        max_len_per_column[j] = max_len;
    }

    for (size_t i = 0; i < grows; ++i)
        for (size_t j = 0; j < gcols; ++j)
            output << (j == 0 ? "\n" : " ") << setw(max_len_per_column[j]) << A.m_matrix[i][j] << (j == gcols - 1 ? "" : "  ");

    cout << '\n';

    return output;
}

istream &operator >> (istream &input, Matrix &A)
{
    cout << endl;
    for (int i = 0; i < A.grows; i++) {
        for (int j = 0; j < A.gcols; j++) {
            input >> A.m_matrix[i][j];
        }
    }
    
    return input;
}

// Returns 3 values
tuple<Matrix, double, int> Matrix::powerIter(int rowNum, double tolerance){
    // Picks a classic X vector
    Matrix X(rowNum,1,1.0);
    // Initiates X vector with values 1,2,3,4
    for (int i = 1; i <= rowNum; i++) {
        X(i-1,0) = i;
    }
    int errorCode = 0;
    double difference = 1.0; // Initiall value greater than tolerance
    int j = 0;
    int location;
    // Defined to find the value between last two eigen values
    vector<double> eigen;
    double eigenvalue = 0.0;
    eigen.push_back(0.0);
    
    while(abs(difference) > tolerance) // breaks out when reached tolerance
    {
        j++;
        // Normalize X vector with infinite norm
        for (int i = 0; i < rowNum; ++i)
        {
            eigenvalue = X(0,0);
            if (abs(X(i,0)) >= abs(eigenvalue))
            {
                // Take the value of the infinite norm as your eigenvalue
                eigenvalue = X(i,0);
                location = i;
            }
        }
        if (j >= 5e5) {
            cout << "Oops, that was a nasty complex number wasn't it?" << endl;
            cout << "ERROR! Returning code black, code black!";
            errorCode = -1;
            return make_tuple(X,0.0,errorCode);
        }
        eigen.push_back(eigenvalue);
        difference = eigen[j] - eigen[j-1];
        // Normalize X vector with its infinite norm
        X = X / eigenvalue;
        
        // Multiply The matrix with X vector
        X = (*this) * X;
    }
    
    // Take the X vector and what you've found is an eigenvector!
    X = X / eigenvalue;
    return make_tuple(X,eigenvalue,errorCode);
}

Matrix Matrix::deflation(Matrix &X, double &eigenvalue)
{
    // Deflation formula exactly applied
    double denominator = eigenvalue / (X.transpose() * X)(0,0);
    Matrix Xtrans = X.transpose();
    Matrix RHS = (X * Xtrans);
    Matrix RHS2 = RHS * denominator;
    Matrix A2 = *this - RHS2;
return A2;
}

Vector::Vector(int rowSize){
    grows = rowSize;
    double init = 0.0;
    m_vector.resize(grows, init);
}

Vector::Vector(int rowSize, double init){
    grows = rowSize;
    m_vector.resize(rowSize, init);
}

// Copy Constructor
Vector::Vector(const Vector &B)
{
    this->grows = B.rows();
    this->m_vector = B.m_vector;
    
}

void Vector::conservativeResize(int rowSize)
{
    Vector B(grows);
    for (int i = 0; i < grows; ++i)
        B(i) = m_vector[i];
    grows = rowSize;
    double init = 0.0;
    m_vector.resize(rowSize, init);

    for (int i = 0; i < grows; ++i)
        this->m_vector[i] = B(i);
}

void Vector::resize(int rowSize)
{
    grows = rowSize;
    double init = 0.0;
    m_vector.resize(rowSize, init);
}

void Vector::random(int rowSize)
{
    grows = rowSize;
    Vect B = Vect::Random(grows);
    double init = 0.0;
    m_vector.resize(rowSize, init);

    for (int i = 0; i < grows; ++i)
        this->m_vector[i] = B(i);
}

Vector::~Vector(){

}

// Addition of Two Matrices
Vector Vector::operator+(Vector &B){
    Vector sum(grows, 0.0);
    int i,j;
    for (i = 0; i < grows; i++)
        sum(i) = this->m_vector[i] + B(i);
    return sum;
}

// Subtraction of Two Matrices
Vector Vector::operator-(Vector & B){
    Vector diff(grows, 0.0);
    int i,j;
    for (i = 0; i < grows; i++)
        diff(i) = this->m_vector[i] - B(i);

    return diff;
}

// Multiplication of Two Matrices
Vector Vector::operator*(Vector & B)
{
    Vector multip(grows,0.0);
    int i;
    double temp = 0.0;
    for (i = 0; i < grows; i++)
        multip(i) = m_vector[i] * B(i);
            
    return multip;
}

// Scalar Addition
Vector Vector::operator+(double scalar){
    Vector result(grows,0.0);
    int i;
    for (i = 0; i < grows; i++)
        result(i) = this->m_vector[i] + scalar;

    return result;
}

// Scalar Subraction
Vector Vector::operator-(double scalar){
    Vector result(grows,0.0);
    int i;
    for (i = 0; i < grows; i++)
        result(i) = this->m_vector[i] - scalar;

    return result;
}

// Scalar Multiplication
Vector Vector::operator*(double scalar){
    Vector result(grows,0.0);
    int i;
    for (i = 0; i < grows; i++)
        result(i) = this->m_vector[i] * scalar;

    return result;
}

// Scalar Division
Vector Vector::operator/(double scalar){
    Vector result(grows,0.0);
    int i;
    for (i = 0; i < grows; i++)
        result(i) = this->m_vector[i] / scalar;

    return result;
}


// Returns value of given location when asked in the form A(x,y)
double& Vector::operator()(const int &rowNo)
{
    return this->m_vector[rowNo];
}

// No brainer - returns row #
int Vector::rows() const
{
    return this->grows;
}

// Prints the Vector beautifully
// ostream& operator << (ostream &output, const Vector &A)
// {
//     cout << endl;
//     for (int i = 0; i < A.grows; i++) 
//             output << A.m_vector[i] << endl;

//     return output;
// }

ostream& operator << (ostream &output, const Vector &A)
{

    int grows = A.grows;
    size_t max_len_per_column[grows];

    for (size_t j = 0; j < grows; ++j) {
        size_t max_len {};

        for (size_t i = 0; i < grows; ++i)
            if (const auto num_length {number_of_digits(A.m_vector[i])}; num_length > max_len)
                max_len = num_length;

        max_len_per_column[j] = max_len;
    }

    for (size_t i = 0; i < grows; ++i)
            output << setw(max_len_per_column[i]) << A.m_vector[i] << endl;

    cout << '\n';

    return output;
}

istream &operator >> (istream &input, Vector &A)
{
    cout << endl;
    for (int i = 0; i < A.grows; i++) 
        input >> A.m_vector[i];
    
    return input;
}