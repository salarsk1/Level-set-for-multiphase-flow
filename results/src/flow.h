#pragma once
#include"calculation.h"
#include"alg.h"
#include<iostream>
#include<fstream>
#include<memory>
using std::vector;
class flow{
public:
	void flow_solver(int, int, double,double, vector<vector<double>>&,vector<vector<double>>&);
};
