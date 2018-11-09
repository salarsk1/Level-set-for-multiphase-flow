#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "vector.h"

#ifdef DEBUGGING_ON
#define DEBUG 1
#else
#define DEBUG 0
#endif

/** Matrix class with row-major format. */
class Matrix {
public:
    //! Constructs an uninitialized matrix with a given size.
    Matrix(int rows, int cols) : _data(rows*cols), _rows(rows), _cols(cols) {}
    //! Constructs a matrix from a nested initializer_list.
    Matrix(const std::initializer_list<std::initializer_list<double>> &x);
    //! Constructs a matrix by copying matrix c. 
    Matrix(const Matrix &c);
    //! Constructs a matrix by moving data from matrix c.
    Matrix(Matrix &&c);

    //! Copies matrix c to this matrix.
    Matrix& operator=(const Matrix &c);
    //! Moves matrix c to this matrix.
    Matrix& operator=(Matrix &&c);

    //! Assigns all elements of this matrix to scalar v.
    void assign(double v);
    
    //! Returns element i,j from this matrix by reference.
    double& operator()(int i, int j);
    //! Returns element i,j from this matrix by value.
    double  operator()(int i, int j) const;

    double *begin()             { return &_data[0]; }
    const double *begin() const { return &_data[0]; }
    double *end()               { return &_data[0] + size(); }
    const double *end() const   { return &_data[0] + size(); }
   
    //! Returns the number of rows in this matrix.
    int rows() const { return _rows; }
    //! Returns the number of columns in this matrix.
    int cols() const { return _cols; }
    //! Returns the number of elements in the matrix. 
    int size() const { return rows()*cols(); }

private:
    //! Stores the data of this matrix in a flat vector. 
    std::vector<double> _data;
    //! Stores the number of rows in the matrix.
    int _rows;
    //! Stores the number of columns in the matrix.
    int _cols;
};
// Constructs a matrix from nested initializer_lists.
inline Matrix::Matrix(const std::initializer_list<
                                std::initializer_list<double>> &x) {
    _rows = x.size();
    if (_rows>0) _cols = x.begin()->size();
    _data.assign(_rows*_cols, 0.0);
    int i=0; 
    for (auto &row: x) {
        // Check that all rows have the same number of columns.
        if (_cols!=row.size()) throw std::string("Invalid matrix");
        int j=0;
        for (auto &v: row) {
            (*this)(i,j++) = v;
        }
        i++;
    }
}
inline Matrix zeros(int rows, int cols) {
    Matrix A(rows, cols);
    A.assign(0.0);
    return A;
}

// Returns a column vector in this matrix.
inline Vector get_col(int cid, const Matrix &A) {
    auto v = zeros(A.rows());
    for (int i=0; i<A.rows(); ++i) {
        v(i) = A(i, cid);
    }
    return v;
}

inline Vector get_row(int rid, const Matrix &A) {
    auto v = zeros(A.cols());
    for (int i=0; i<A.cols(); ++i) {
        v(i) = A(rid, i);
    }
    return v;
}

// Assign a vector to a matrix column. 
inline void assign_col(int cid, const Vector &a, Matrix &A) {
    if (DEBUG && (A.cols() <= cid || a.size() != A.rows())) {
        throw std::string("Bad matrix column assignment.");
    }
    for (int i=0; i<a.size(); ++i) {
        A(i, cid) = a(i);
    }
}

// Assign a vector to a matrix row. 
inline void assign_row(int rid, const Vector &a, Matrix &A) {
    if (DEBUG && (A.rows() <= rid || a.size() != A.cols())) {
        throw std::string("Bad matrix row assignment.");
    }
    for (int i=0; i<a.size(); ++i) {
        A(rid, i) = a(i);
    }
}

// Returns a matrix with random value between 0 and 1.
inline Matrix rand(int rows, int cols) {
    Matrix A(rows, cols);
    for (int i=0; i<rows; ++i)
        for (int j=0; j<cols; ++j)
            A(i,j) = rand_double();
    return A;
}

// Returns an identity matrix.
inline Matrix eye(int rows) {
    Matrix A = zeros(rows, rows);
    for (int i=0; i<rows; ++i)
        A(i,i) = 1.0;
    return A;
}

inline Matrix::Matrix(const Matrix &c) {
    _data = c._data;
    _rows = c._rows;
    _cols = c._cols;
}

inline Matrix::Matrix(Matrix &&c) {
    _data = std::move(c._data);
    _rows = std::move(c._rows);
    _cols = std::move(c._cols);
}

inline Matrix& Matrix::operator=(const Matrix &c) {
    _data = c._data;
    _rows = c._rows;
    _cols = c._cols;
    return *this;
}

inline Matrix& Matrix::operator=(Matrix &&c) {
    _data = std::move(c._data);
    _rows = std::move(c._rows);
    _cols = std::move(c._cols);
    return *this;
}

inline void Matrix::assign(double v) {
    for (auto x=_data.begin(); x!=_data.end(); ++x) *x = v;
}

inline double& Matrix::operator()(int i, int j) {
    if (DEBUG && (i>=_rows || j>=_cols)) {
        throw std::string("Matrix indexing out of bounds.");
    }
    return _data[j+i*_cols]; 
}
inline double Matrix::operator()(int i, int j) const { 
    if (DEBUG && (i>=_rows || j>=_cols)) {
        throw std::string("Matrix indexing out of bounds.");
    }
    return _data[j+i*_cols]; 
}

inline Matrix& operator+=(Matrix &A, const Matrix &B) {
    if(DEBUG && A.cols()!=B.cols())
        throw std::string("Matrix-Matrix add - columns mismatch.");
    if(DEBUG && A.rows()!=B.rows())
        throw std::string("Matrix-Matrix add - rows mismatch.");
    for (int i=0; i<A.rows(); ++i)
        for (int j=0; j<A.cols(); ++j)
            A(i,j) += B(i,j); 
    return A;
}

inline Matrix operator+(const Matrix &A, const Matrix &B) {
    Matrix C(A);
    C+=B;
    return C;
}

inline Matrix operator+(const Matrix &A, const double &s) {
    Matrix C(A);
    for (auto &c: C) c += s;
    return C;
}

inline Matrix& operator-=(Matrix &A, const Matrix &B) {
    if(DEBUG && A.cols()!=B.cols())
        throw std::string("Matrix-Matrix subtract - columns mismatch.");
    if(DEBUG && A.rows()!=B.rows())
        throw std::string("Matrix-Matrix subtract - rows mismatch.");
    for (int i=0; i<A.rows(); ++i)
        for (int j=0; j<A.cols(); ++j)
            A(i,j) -= B(i,j); 
    return A;
}

inline Matrix operator-(const Matrix &A, const Matrix &B) {
    Matrix C(A);
    C -= B;
    return C;
}

inline Matrix operator-(const Matrix &A, const double &s) {
    Matrix C(A);
    for (auto &c: C)
        c -= s;
    return C;
}

//! Returns a vector, vector cross-product.
inline Matrix operator*(const Vector &a, const Vector &b) {
	Matrix C(a.size(), b.size());
	for (int i=0;i<a.size();i++)
		for (int j=0;j<b.size();j++)
			C(i,j) = a(i)*b(j);
	return C;
}

//! Returns a Matrix, scalar product.
inline Matrix operator*(const Matrix &A, const double b) {
	Matrix C(A.rows(), A.cols());
	for (int i=0;i<C.rows();i++)
		for (int j=0;j<C.cols();j++)
			C(i,j) = A(i,j)*b;
	return C;
}

inline Vector operator*(const Matrix &A, const Vector &x) {
    if(DEBUG && A.cols()!=x.size()){
        std::cout << "Matrix-Vector multiply - bad sizes.\n";
        std::cout << "A is " << A.rows() << " by " << A.cols() << ".\n";
        std::cout << "x is " << x.size() << ".\n";
        exit(1);
    }
    Vector y(A.rows());
    for (int i=0; i<A.rows(); ++i) {
        y(i) = A(i,0)*x(0);
        for (int j=1; j<A.cols(); ++j) 
            y(i) += A(i,j) * x(j); 
    }
    return y;
}

//! Returns a matrix, vector product.
inline Vector operator*(const Vector &x, const Matrix &A) {
    if(DEBUG && x.size()!=A.rows()){
        std::cout << "Vector-Matrix multiply - bad sizes.\n";
        std::cout << "x is " << x.size() << ".\n";
        std::cout << "A is " << A.rows() << " by " << A.cols() << ".\n";
        throw std::string("");
    }
    Vector y(A.cols());
    for (int j=0; j<A.cols(); ++j) {
        y(j) = x(0)*A(0,j);
        for (int i=1; i<A.rows(); ++i)
            y(j) += x(i)*A(i,j);
    }
    return y;
}

//! Calls LAPACK lu-decomposition routines.
extern "C" {
void dgetrf_(int*, int*, double*, int*, int*, int*);
void dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
}
inline Vector solve(const Matrix &A, const Vector &y) {
    char t = 'N';
    int dim[] = {A.rows()};
    int nrhs[] = {1};
    int *LDA=dim, *LDB=dim;
    int info;

    std::vector<int> ipiv(*dim);
    Matrix LU(A);
    Vector x(y);
    dgetrf_(dim, dim, LU.begin(), LDA, ipiv.data(), &info);
    dgetrs_(&t, dim, nrhs, LU.begin(), LDA, ipiv.data(), x.begin(), LDB, &info);
    return x;
}

//! Returns a matrix product.
inline Matrix operator*(const Matrix &A, const Matrix &B) {

    if (DEBUG && A.cols()!=B.rows()) {
        std::cout << "Matrix multiply - bad sizes.\n";
        std::cout << "A is " << A.rows() << " by " << A.cols() << ".\n";
        std::cout << "B is " << B.rows() << " by " << B.cols() << ".\n";
        throw std::string("");
    }

    Matrix C(A.rows(), B.cols());
    for (int i=0; i<A.rows(); ++i) {
        for (int j=0; j<B.cols(); ++j) {
            C(i,j) = A(i,0) * B(0,j);
            for (int k=1; k<A.cols(); ++k) {
                C(i,j) += A(i,k) * B(k, j);
            }
        }
    }
    return C;
}

//! Returns the determinate of a matrix.
inline double det(const Matrix &A) {
    if (A.rows() == 3 && A.cols() == 3) {
        return A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1)) 
             + A(0,1)*(A(1,2)*A(2,0)-A(1,0)*A(2,2))
    		 + A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
    }
    else if (A.rows() == 2 && A.cols() == 2) {
        return A(0,0)*A(1,1) - A(0,1)*A(1,0);
    }
    else std::cout << "matrix det - not 3x3 or 2x2\n";
    return 0;
}

//! Returns the trace of a matrix.
inline double trace(const Matrix &A) {
    if (DEBUG && A.cols()!=A.rows()) {
        throw std::string("Matrix trace: not square matrix.");
    }
    double sum = 0.0;
    for (int i=0; i<A.rows(); ++i) {
        sum += A(i,i);
    }
    return sum;
}

//! Returns the sum of all elements (or absolute values) in a matrix.
inline double sum(const Matrix &A, bool is_abs_sum=false) {
    double s=0.0;
    for (int i=0; i<A.rows(); ++i) {
        for (int j=0; j<A.cols(); ++j) {
            if (is_abs_sum) s += fabs(A(i,j));
            else s += A(i,j);
        }
    }
    return s;
}

//! Returns the inverse of a matrix.
inline Matrix inv(const Matrix &A) {
    double invDetA = 1.0 / det(A);
    if (A.rows() == 3 && A.cols() == 3) {
        Matrix Ai(3,3);
		Ai(0,0) = (A(1,1)*A(2,2)-A(1,2)*A(2,1))*invDetA;
		Ai(0,1) = (A(0,2)*A(2,1)-A(0,1)*A(2,2))*invDetA;
		Ai(0,2) = (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invDetA;
		Ai(1,0) = (A(1,2)*A(2,0)-A(1,0)*A(2,2))*invDetA;
		Ai(1,1) = (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invDetA;
		Ai(1,2) = (A(0,2)*A(1,0)-A(0,0)*A(1,2))*invDetA;
		Ai(2,0) = (A(1,0)*A(2,1)-A(1,1)*A(2,0))*invDetA;
		Ai(2,1) = (A(0,1)*A(2,0)-A(0,0)*A(2,1))*invDetA;
		Ai(2,2) = (A(0,0)*A(1,1)-A(0,1)*A(1,0))*invDetA;
        return Ai;
    }
    else if (A.rows() == 2 && A.cols() == 2) {
        Matrix Ai(2,2);
        Ai(0,0) =  A(1,1) * invDetA;
        Ai(0,1) = -A(0,1) * invDetA;
        Ai(1,0) = -A(1,0) * invDetA;
        Ai(1,1) =  A(0,0) * invDetA;
        return Ai;
    }
    else std::cout << "matrix inverse - not 3x3 or 2x2\n";
    return A;
}

//! Returns the transpose of a matrix.
inline double norm(const Matrix &A) {
    double r = 0.0;
    for (auto x: A) r += x*x;
    return sqrt(r);
}

//! Returns the transpose of a matrix. 
inline Matrix transpose(const Matrix &A) {
    Matrix B(A.cols(), A.rows());
    for (int i=0; i<A.rows(); ++i)
        for (int j=0; j<A.cols(); ++j)
            B(j,i) = A(i,j);
    return B;
}

//! Outputs a Matrix to an output stream like cout.
static std::ostream& operator<<(std::ostream &o, const Matrix &A) {
    for (int i=0; i<A.rows(); ++i) {
        for (int j=0; j<A.cols(); ++j) {
            if (j) o << " ";
            o << A(i,j);
        }
        o << "\n";
    }
    return o;
}

