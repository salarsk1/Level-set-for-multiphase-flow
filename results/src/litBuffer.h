//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include "parallel.h"
class litBuffer{
public:
	void litBuffer_m_init();
	bool*    getL1Buffer(size_t &n, char band='F', int* i2b = nullptr);
	int*     getI1Buffer(size_t &n, char band='F', int* i2b = nullptr);
	double*  getR1Buffer(size_t &n, char band='F', int* i2b = nullptr);
	double** getR2Buffer(size_t&, size_t&, char band='F', int nr2=0, int* i2b=nullptr);
	void freeL1Buffer(bool    *l1);
	void freeI1Buffer(int     *i1);
	void freeR1Buffer(double  *r1);
	void freeR2Buffer(size_t, size_t, double **r2);
};
