//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include"parallel.h"
#include<string>
#include<vector>
#include<algorithm>
#include<fstream>
#include"string.h"
using std::string;
using std::vector;
using std::fstream;
class param_t;
class litparam{
public:
	void init_param(int argc, char *argv[], const string& input_filrname="\0");
	void destroy_param();
	void destroy_param_recursive(param_t *param);
	void parse_commandline_param(int argc, char * argv[], int &ierr);
	param_t* get_new_param();
	void build_mpi_filename(string &mpi_filename,const string filename);
	int get_next_token(int,int,int&,int&,int&,int&,int&,int&,int&,int&,string&,bool&,string,\
							 MPI_File, MPI_Status, int);
	void read_param_file(std::string filename);
	void add_param_line(const string line, int argc, int* argv);
	param_t* get_param(const string name);
	void set_next_param(param_t* param, const string name);
	int get_param_type(const param_t &param, const int iarg);
	void dump_param();
	void dump_param_usage();
	bool check_param(const string name);
	bool get_logical_param(const string name, bool defa, int arg=1);
	int get_integer_param(const string name, int defa, int arg=1);
	double get_real_param(const string name, double defa, int arg = 1);
	string get_string_param(const string name, string defa, int arg);
	std::array<int,3> get_integer3_param(const string name, std::vector<int> defa, int arg=1);
	std::array<double,3>get_real3_param(const string name,const std::vector<double> defa,int arg=1);
	void set_param(string &value, const string name, int index=-1,const string defa="");
	void set_param_token(string &value,const param_t &param, int index=-1);
	inline std::string trim(std::string s, std::string ws=" \t\n\r") {
		if (s.empty()) return s;
		size_t a=s.find_first_not_of(ws);
		size_t b = s.find_last_not_of(ws)+1;
		return s.substr(a, b-a);
	}
	// split the string like a python function
	inline std::vector<std::string> split(std::string s, std::string delims=" \t\n\r"){
		std::vector<std::string> pieces;
	    size_t begin=0, end=0;
	    while (end != std::string::npos) {
	        begin = s.find_first_not_of(delims, end);
	        end   = s.find_first_of(delims, begin);
	        if (begin != std::string::npos)
	            pieces.push_back(trim(s.substr(begin, end-begin)));
	    }
	    return pieces;
	}
	inline std::string read_line(std::fstream &fid) {
		std::string line;
		std::getline(fid, line);
		return line;
	}

	// returns the lowercase string
	inline string lower(string s) {
		std::transform(s.begin(), s.end(), s.begin(),
		static_cast<int(*)(int)> (tolower));
		return s;
	}	
private:
	const int NUMERICAL_PARAM_TYPE=1;
	const int LOGICAL_PARAM_TYPE=2;
	const int STRING_PARAM_TYPE=3;
	const int NO_PARAM_TYPE=4;
	const int MAX_PARAM_LEN=64;
	const int MAX_STRING_LEN=256;
	const int MAX_LINE_LEN=256;
	const int MAX_RECURSION_LEVEL=8;
	const int MAX_TOKEN_LEN=64;
	const int MAX_TOKEN_COUNT=32;
	int recursion_level = 0;
	static param_t* first_param_ptr;
};
class param_t{
public:
	std::string line;
	int argc;
	int *argv;
	int request_count;
	param_t *next;
};

