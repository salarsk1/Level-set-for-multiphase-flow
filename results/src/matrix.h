//Written by Salar Safarkhani

#include<vector>
#include<iostream>
#include<algorithm>
#include"vector.h"
#include<initializer_list>
class Matrix{
public:
Matrix(int rows, int cols) :_data(rows*cols), _rows(rows), _cols(cols){}
Matrix(const std::initializer_list<std::initializer_list<double>> &x);
Matrix(const Matrix &c);
Matrix(Matrix &&c);
Matrix &operator=(const Matrix &c);
Matrix &operator=(Matrix &&c);
double &operator()(int i, int j);
double  operator()(int i, int j) const;
double *begin()             {return &_data[0];}
const double *begin() const {return &_data[0];}
double *end()               {return &_data[0] + size();}
const double *end() const   {return &_data[0] + size();}
int rows() const {return _rows;}
int cols() const {return _cols;}
int size() const {return rows()*cols();}
void assign(double v);
private:
std::vector<double> _data;
int _rows;
int _cols;
};
inline Matrix::Matrix(const std::initializer_list<\
  									 std::initializer_list<double>> &x){
_rows = x.size();
if(_rows>0) _cols = x.begin()->size();
_data.assign(_rows*_cols, 0.0);
int i=0;
for(auto &row:x){
	int j=0;
	for(auto &v:row){
		(*this)(i, j++) = v;
	}
	i++;
};
};
inline Matrix::Matrix(const Matrix &c){
_data = c._data;
_rows = c._rows;
_cols = c._cols;
};
inline Matrix::Matrix(Matrix &&c){
_data = c._data;
_rows = c._rows;
_cols = c._cols;
};
inline Matrix &Matrix::operator=(const Matrix &c){
_data = c._data;
_rows = c._rows;
_cols = c._cols;
return *this;
};
inline Matrix &Matrix::operator=(Matrix &&c){
_data = std::move(c._data);
_rows = std::move(c._rows);
_cols = std::move(c._cols);
return *this;
};

inline double& Matrix::operator()(int i, int j){
    if (i>=_rows || j>=_cols) {
        throw std::string("Matrix indexing out of bounds.");
    }
    return _data[j+i*_cols];
}
inline double Matrix::operator()(int i, int j) const{
    if (i>=_rows || j>=_cols) {
        throw std::string("Matrix indexing out of bounds.");
    }
    return _data[j+i*_cols];
}

inline Matrix& operator+=(Matrix &A, const Matrix &B) {
    if(A.cols()!=B.cols())
        throw std::string("Matrix-Matrix add - columns mismatch.");
    if(A.rows()!=B.rows())
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
    if(A.cols()!=B.cols())
        throw std::string("Matrix-Matrix subtract - columns mismatch.");
    if(A.rows()!=B.rows())
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

//! Returns a Matrix, scalar product.
inline Matrix operator*(const Matrix &A, const double b) {
   Matrix C(A.rows(), A.cols());
   for (int i=0;i<C.rows();i++)
      for (int j=0;j<C.cols();j++)
         C(i,j) = A(i,j)*b;
   return C;
}

inline Vector operator*(const Matrix &A, const Vector &x) {
    if(A.cols()!=x.size()){
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

inline Vector operator*(const Vector &x, const Matrix &A) {
    if(x.size()!=A.rows()){
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

inline Matrix operator*(const Matrix &A, const Matrix &B) {

    if (A.cols()!=B.rows()) {
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

