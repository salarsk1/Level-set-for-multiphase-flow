//Written by Salar Safarkhani

#pragma once
#include "datag.h"
#include"parallel.h"
#include"param.h"
class lit_mfile_t;
class monitor{
public:
	void lit_monitor_set_onevalue(const int i, const double v);
	void lit_monitor_dump_values_step();
	void lit_monitor_dump_values_iter(const int ntype, const int iter);
	void lit_monitor_set_single_values(const double v1=0,const double v2=0,const double v3=0,\
		  const double v4=0,const double v5=0,const double v6=0,const double v7=0,const double v8=0,\
		  const double v9=0,const double v10=0,const double v11=0,const double v12=0,\
		  const double v13=0,const double v14=0,const double v15=0,const double v16=0,\
		  const double v17=0,const double v18=0,const double v19=0);
	void lit_monitor_pre_init();
	void lit_monitor_post_init();
	void lit_monitor_create_gnuplot(const int iy1=0,const int iy2=0,\
											  const int iy3=0,const int iy4=0,\
											  const int iy5=0,const int iy6=0,\
											  const int iy7=0,const int iy8=0,\
											  const int iy9=0,const bool send=false);
	void lit_monitor_create_file_step(const std::string filename, const int ncols);
	void lit_monitor_create_file_iter(const std::string filename,const int ncols,\
			  								    const int ntype);
	void lit_monitor_select_file(const std::string filename);
	void lit_monitor_set_header(const int icol, const std::string header, \
										 const char col_type);
	void lit_monitor_set_array_values(double *val);
	void lit_monitor_log(const std::string text);
	void lit_monitor_finalize();

	static double monitor_gradG_max, monitor_gradG_min;
	static int monitor_reinit_steps, monitor_nSubCycles;
private:
	const int lit_nfiles_max = 32;
	const int col_len = 14; //need not to be defined; checked
	static int lit_nfiles;
	static int lit_ifile;
	int lit_file_log;

	static lit_mfile_t *lit_mfiles;

	static std::string gfile;
};
class lit_mfile_t{
public:
	int iunit;
	std::string filename;
	int freq;
	int ncols;
	std::string *header;
	std::string col_type;
	//char *col_type; // checked
	double *val;
};
