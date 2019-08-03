#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <memory>
#include <cmath>

#ifdef DEBUGGING_ON
#define DEBUG 1
#else
#define DEBUG 0
#endif

/** A simple vector class for linear algebra operations. */
class Vector {
public:
    //! Constructs an uninitialized vector of length, len.
    Vector(int len=0);
    //! Constructs a vector from an initializer list.
    Vector(std::initializer_list<double> x);
    //! Constructs by copying another vector.
    Vector(const Vector &r);
    //! Constructs by moving data from another vector.
    Vector(Vector &&r);

    //! Assigns a vector so that all elements are a constant value.
    void assign(double v);

    //! Returns the ith element in the vector by reference. 
    double &operator()(int i)       { return _members[i]; }
    //! Returns the ith element in the vector by value. 
    double  operator()(int i) const { return _members[i]; }
    //! Returns a subset of elements of this vector.
    Vector operator()(const std::vector<int> &set) const;

    //! Returns the count of elements in the vector.
    int size() const { return _size; }
    bool empty() const { return _size == 0; }

    //! Returns pointer to the first element.
    double *begin() { return _members.get(); }
    //! Returns pointer to the end element.
    double *end()   { return _members.get() + _size; }
    //! Returns const pointer to the first element.
    double const* begin() const { return _members.get(); }
    //! Returns const pointer to the end element.
    double const* end()   const { return _members.get() + _size; }

    //! Moves the data from vector r to this vector, leaving r empty.
    Vector& operator=(Vector &&r);
    //! Assigns the vector to another vector, r.
    Vector& operator=(const Vector &r);
    //! Scales the vector by a value.
    Vector& operator*=(double x);
    //! Divide the vector by a value.
    Vector& operator/=(double x);

    //! Write the vector to a ascii file (one element per line).
    void write_ascii(const char *path) const;

private:
    //! Pointer to first element in the vector.
    std::unique_ptr<double[]> _members;
    //! Number of elements in the vector.
    int     _size;
};

//! Returns a vector initialized to all zeros.
inline Vector zeros(int len) {
    Vector x(len);
    x.assign(0.0);
    return x;
}

inline Vector::Vector(int len)
  : _members(new double[len]), _size(len) {}

inline Vector::Vector(std::initializer_list<double> x) 
  : _members(new double[x.size()]), _size(x.size()) {
    std::copy(x.begin(), x.end(), _members.get());
}

// Moves c's data to here.
inline Vector::Vector(Vector &&c) {
    _size = c.size();
    _members = std::move(c._members);
}

// Makes a copy of c.
inline Vector::Vector(const Vector &c) {
    _size = c.size();
    _members.reset(new double[c.size()]);
    std::copy(c.begin(), c.end(), _members.get());
}

inline Vector& Vector::operator=(const Vector &r) {
    if (_size != r.size()) {
        _size = r.size();
        _members.reset(new double[r.size()]);
    }
    std::copy(r.begin(), r.end(), _members.get());
    return *this;
}

//! Moves data from r to this vector.
inline Vector& Vector::operator=(Vector &&r) {
    _members = std::move(r._members);
    _size    = r.size();
    return *this;
}

inline Vector& Vector::operator*=(double x) {
    for (auto p=begin(); p!=end(); ++p) *p *= x;
    return *this;
}

inline Vector& Vector::operator/=(double x) {
    for (auto p=begin(); p!=end(); ++p) *p /= x;
    return *this;
}

inline void Vector::assign(double v) {
    for (auto p=begin(); p!=end(); ++p) *p = v;
}

inline Vector Vector::operator()(const std::vector<int> &set) const {
    Vector r(set.size());
    for (size_t i=0; i<set.size(); ++i) r(i) = (*this)(set[i]);
    return r;
}

//! Outputs the vector to an output stream like cout.
static std::ostream& operator<<(std::ostream &o, const Vector &x) {
    for (int i=0; i<x.size(); ++i) {
        if (i) o << "   ";
        o << std::setw(12) << x(i);
    }
    return o << "\n";
}

//! Outputs the vector to an output file.
inline void Vector::write_ascii(const char *path) const {
    std::fstream fid(path, std::ios::out);
    for (int row=0; row<size(); ++row) {
        fid << std::setprecision(16)<< (*this)(row) << "\n";
    }
} 

//! Returns the dot product of vectors a and b.
inline double dot(const Vector &x, const Vector &y) {
    if (DEBUG && x.size()!=y.size()) 
        throw std::string("Vector dot product: size mismatch");
    double result = 0.0;
    for (int i=0; i<x.size(); ++i) {
        result += x(i)*y(i);
    }
    return result;
}

//! Returns the cross product of 3D vectors a and b.
inline Vector cross(const Vector &a, const Vector &b) {
    if (DEBUG && (a.size()!=b.size() || a.size()!=3)) 
        throw std::string("Vector cross product: size is not 3");
    Vector c(3);
    c(0) = a(1)*b(2) - a(2)*b(1);
    c(1) = a(2)*b(0) - a(0)*b(2);
    c(2) = a(0)*b(1) - a(1)*b(0);
    return c;
}

//! Returns the norm of a vector v.
template <typename T>
inline double norm(const T &v) {
    double result = 0.0;
    for (auto x: v) result += x*x;
    if (DEBUG && result != result)  {
        throw std::string("Vector norm is NaN");
    }        
    return sqrt(result);
}

//! Returns the max element of vector v.
inline double max(const Vector &v) {
    if (DEBUG && v.empty()) 
        throw std::string("Can't take max of empty vector");
    return *std::max_element(v.begin(), v.end());
}

inline double min(const Vector &v) {
    if (DEBUG && v.empty()) 
        throw std::string("Can't take min of empty vector");
    return *std::min_element(v.begin(), v.end());
}

//! Returns the max element position of vector v.
inline int max_loc(const Vector &v) {
    if (DEBUG && v.empty()) 
        throw std::string("Can't take max of empty vector");
    return std::max_element(v.begin(), v.end()) - v.begin();
}

inline int min_loc(const Vector &v) {
    if (DEBUG && v.empty()) 
        throw std::string("Can't take min of empty vector");
    return std::min_element(v.begin(), v.end()) - v.begin();
}

//! Keeping track of original indexes for sorted vector. 
inline std::vector<int> sort_indexes(const Vector &v) {
    std::vector<int> index(v.size());
    // Initialize original index locations.
    for (int i=0; i<v.size(); ++i) index[i] = i;

    // Sort indexes based on comparing values in the given vector. 
    std::sort(index.begin(), index.end(), 
              [&v](int a, int b) {return v(a)<v(b);});

    return index;
}

//! Returns the sum of vector v.
inline double sum(const Vector &v) {
    double result = 0.0;
    for (auto x: v) result += x;
    return result;
}

//! Returns the absolute of vector v.
inline Vector abs(const Vector &v) {
	Vector result(v.size());
	for (int i=0;i<v.size();i++) result(i) = fabs(v(i));
	return result;
}

//! Returns sign function result.
inline double sign(const double s) {
	if (s>0) return 1;
	else if (s<0) return -1;
	else return 0;
}

//! Returns sign function of a vector.
inline Vector sign(const Vector &v) {
    Vector sv(v);
    for (auto &s : sv) {
        if (s>0)      s = 1.0;
        else if (s<0) s = -1.0;
        else          s = 0.0;
    }
    return sv;
}

//! Adds right hand side vector into the left hand side.
inline Vector& operator+=(Vector &a, const Vector &b) {
    if (DEBUG && a.size()!=b.size()) 
        throw std::string("Vector addition: size mismatch");
    for (int i=0; i<a.size(); ++i) a(i) += b(i);
    return a;
}

//! Returns the sum of vectors x and y.
inline Vector operator+(Vector x, const Vector &y) {
    return x += y;
}

//! Subtracts right hand side vector from the left hand side.
inline Vector& operator-=(Vector &a, const Vector &b) {
    if (DEBUG && a.size()!=b.size()) 
        throw std::string("Vector subtraction: size mismatch");
    for (int i=0; i<a.size(); ++i) a(i) -= b(i);
    return a;
}

//! Returns the difference of vectors x and y.
inline Vector operator-(Vector x, const Vector &y) {
    return x -= y;
}

//! Adds a scalar y to all elements of x and returns a new vector.
inline Vector operator+(const Vector &x, double y) {
    Vector c(x);
    for (int i=0; i<c.size(); ++i) c(i)+=y;
    return c;
}
inline Vector operator+(double y, const Vector &x) {
    return x + y;
}

//! Subtracts a scalar y to all elements of x and returns a new vector.
inline Vector operator-(const Vector &x, double y) {
    Vector c(x);
    for (int i=0; i<c.size(); ++i) c(i)-=y;
    return c;
}
inline Vector operator-(double y, const Vector &x) {
    return x - y;
}

//! Returns a copy of x scaled by s.
inline Vector operator*(const Vector &x, double s) {
    Vector c(x);
    for (int i=0; i<c.size(); ++i) c(i)*=s;
    return c;
}
inline Vector operator*(double s, const Vector &x) {
    return x * s;
}

//! Creates vector of all ones. 
inline Vector ones(int size) {
    Vector c(size);
    for (auto &x: c) x = 1.0;
    return c;
}

//! Returns a vector of num evenly spaced points between start and stop.
inline Vector linspace(double start, double stop, int num) {
    double dt = (stop-start)/double(num-1);
    Vector x(num);
    x(0)     = start;
    x(num-1) = stop;
    for (auto ptr=x.begin()+1; ptr!=x.end()-1; ++ptr) {
        *ptr = *(ptr-1) + dt;
    }
    return x;
}

//! Returns a random number between 0 and 1.
inline double rand_double() { return double(rand()) / double(RAND_MAX); }
//! Returns a vector full of random numbers between 0 and 1 (like MATLAB).
inline Vector rand(int n) {
    Vector r(n);
    for (auto &x: r) x = rand_double();
    return r;
}


