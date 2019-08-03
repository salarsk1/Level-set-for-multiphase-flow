//Written by Salar Safarkhani

#pragma once
#include<fstream>
#include<sstream>
#include "datag.h"
#include "timing.h"
#include "parallel.h"
#include "param.h"
#include "monitor.h"
#include "bound.h"
#include "toolbox.h"
#include "bl.h"
#include "litBuffer.h"
#include "reinit.h"
#include "sg.h"
#include "band.h"
#include "advection.h"
#include "gnodes.h"
class io{
public:
	void io_init();
	void read_input(int argc, char *argv[], std::string input_filename, int arg = 0);
	void dumpSolution(bool newv, bool force, int arg = 0);
	void dumpEnsight(const std::string name, const std::string band, bool newv, int arg=2);
	void writeEnsightHeaderScalar(const std::string varname,std::string &,const std::string, \
									MPI_File &, std::string &,std::string&,int &,bool&, const std::string);
	void writeEnsightHeaderVector(const std::string varname,std::string &,const std::string, \
									MPI_File &, std::string &,std::string&,int &,bool&, const std::string);
	void writeEnsightScalar(size_t ,double* var, MPI_File &, int &);
	void writeEnsightVector(size_t, double *var, int id, MPI_File &, int &, int &);
	void dumpEnsightCell(const std::string, const std::string band, const bool newv, int arg =3);
	void writeCaseFileHeader(const std::string name);
	inline bool isVector(const std::string varName){
		std::unique_ptr<litparam> ppar(new litparam);
		if(ppar->trim(varName)=="V" || ppar->trim(varName)=="DIJK") return true;
		else return false;
	};
	
	inline std::vector<string> split(string str){
	    string buf;
	    std::stringstream ss(str);
	    vector<string> tokens;
	    while (ss >> buf) tokens.push_back(buf);
		 return tokens;
	}

	void writeCaseFileCellHeader(const std::string name);
	void write_new_case_file(const std::string name, const int typ);
	void append_case_file(std::string name, const int typ);
	void dump_restart(const bool keep, int arg = 0);
	inline bool exists_file(const std::string &name){
		std::ifstream f(name);
		if(f.good()){
			f.close();
			return true;
		}
		else{
			f.close();
			return false;
		}
	}
	static int dump_nfiles;
	static int nVarDump;
	static int nVarCellDump;
	static std::vector<std::string> varDumpName;
	static std::vector<std::string> varCellDumpName;

private:
	inline double round(double d){
		return floor(0.5+d);
	}
	static string global_file;
};

