//Written by Salar Safarkhani

#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>
#include <cmath>
class Vector{
public:
Vector(int len=0);
Vector(std::initializer_list<double> x);
Vector(const Vector &r);
Vector(Vector &&r);
double &operator()(int i) {return _members[i]; };
double  operator()(int i) const {return _members[i]; };
int size() const {return _size; };
Vector &operator=(Vector &&);
Vector &operator=(const Vector &r);
int size() {return _size;};
double *begin() {return _members.get(); };
double *end()   {return _members.get() + _size; };
double const *begin() const {return _members.get(); };
double const *end()   const {return _members.get() + _size; };
private:
std::unique_ptr<double[]> _members;
int _size;
};
inline Vector::Vector(int len):
	_members(new double[len]), _size(len) {};

inline Vector::Vector(Vector &&c) {
_size = c.size();
_members = std::move(c._members);
};

inline Vector::Vector(const Vector &c) {
_size = c.size();
//_members.reset(new double[c.size()]);
std::copy(c.begin(), c.end(), _members.get());
};
	
inline Vector::Vector(std::initializer_list<double> x):
	_members(new double[x.size()]), _size(x.size()){
	std::copy(x.begin(), x.end(), _members.get());
};

inline Vector& Vector::operator=(const Vector &r) {
    if (_size != r.size()) {
        _size = r.size();
        _members.reset(new double[r.size()]);
    }
    std::copy(r.begin(), r.end(), _members.get());
    return *this;
}

inline Vector& Vector::operator=(Vector &&r) {
    _members = std::move(r._members);
    _size    = r.size();
    return *this;
}






