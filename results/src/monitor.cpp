//Written by Salar Safarkhani

#include "monitor.h"
#include "parallel.h"
#include<iostream>
#include<fstream>
#include<ctime>
#include<boost/filesystem.hpp>
void monitor::lit_monitor_set_onevalue(const int i, const double v){
std::unique_ptr<parallel>pp(new parallel);
if(i>lit_mfiles[lit_ifile-1].ncols){
	std::cout<<"fname = "<<lit_mfiles[lit_ifile-1].filename<<std::endl;
	std::cout<<"i,ncols = "<<i<<lit_mfiles[lit_ifile-1].ncols<<std::endl;
	pp->litError("monitor_set", "too many values");
}
else{
	lit_mfiles[lit_ifile-1].val[i-1] = v;
};
};

void monitor::lit_monitor_dump_values_step(){
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<litparam>ppar(new litparam);
std::string line(4096,' ');
std::string strtemp;
std::ofstream myfile;
myfile.open(gfile);
std::string S1;
if(pp->myrank != 0) return;
for(lit_ifile = 0; lit_ifile<lit_nfiles; lit_ifile++){
	if(lit_mfiles[lit_ifile].freq == 1){
		line.replace(0,  std::to_string(pgv->fs_step).size(),std::to_string(pgv->fs_step));//format
		line.replace(1*col_len,std::to_string(pgv->time).size(),std::to_string(pgv->time));//format
		int offset = 2*col_len;
		for(auto i=0; i<lit_mfiles[lit_ifile].ncols; ++i){
			switch(lit_mfiles[lit_ifile].col_type[i]){
				case 'i':
					S1 = std::to_string(static_cast<int>(lit_mfiles[lit_ifile].val[i]));
					line.replace(offset+i*col_len,S1.size(), S1); //check format
					break;
				case 'r':
					S1 = std::to_string(lit_mfiles[lit_ifile].val[i]);
					line.replace(offset+i*col_len, S1.size(), S1);//check format
					break;
			};
		};
		myfile<<ppar->trim(line); // check this
	}
}
};
void monitor::lit_monitor_dump_values_iter(const int ntype, const int iter){
std::string line(4096, ' ');
std::string strtemp;
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<litparam>ppar(new litparam);
std::ofstream myfile;
myfile.open(gfile);
if(pp->myrank != 0) return;
for(auto i=0; i<lit_nfiles; i++){
	if(lit_mfiles[i].freq == 1){
		line.insert(0, std::to_string(pgv->fs_step)); // check format
		line.insert(1*col_len, std::to_string(pgv->time)); // check format
		line.insert(2*col_len, std::to_string(iter)); // check format
		auto offset = 3*col_len;
		for(auto j=0; j<lit_mfiles[i].ncols; j++){
			if(lit_mfiles[i].col_type[j] == 'i'){
				line.insert(offset+1+j, std::to_string(static_cast<int>(lit_mfiles[i].val[j]))); // format
			}
			else if(lit_mfiles[i].col_type[j] == 'r'){
		line.insert(offset+1+j, std::to_string(lit_mfiles[i].val[j])); // check format
			}
			myfile<<ppar->trim(line);	
		}
	}
}
}
void monitor::lit_monitor_set_single_values(const double v1,const double v2,const double v3,\
		const double v4,const double v5,const double v6,const double v7,const double v8,\
		const double v9,const double v10,const double v11,const double v12,const double v13,\
		const double v14,const double v15,const double v16,const double v17,const double v18,\
		const double v19){
/*if(v1 != 0.0)  lit_monitor_set_onevalue(1 , v1);
if(v2 != 0.0)  lit_monitor_set_onevalue(2 , v2);
if(v3 != 0.0)  lit_monitor_set_onevalue(3 , v3);
if(v4 != 0.0)  lit_monitor_set_onevalue(4 , v4);
if(v5 != 0.0)  lit_monitor_set_onevalue(5 , v5);
if(v6 != 0.0)  lit_monitor_set_onevalue(6 , v6);
if(v7 != 0.0)  lit_monitor_set_onevalue(7 , v7);
if(v8 != 0.0)  lit_monitor_set_onevalue(8 , v8);
if(v9 != 0.0)  lit_monitor_set_onevalue(9 , v9);
if(v10 != 0.0) lit_monitor_set_onevalue(10, v10);
if(v11 != 0.0) lit_monitor_set_onevalue(11, v11);
if(v12 != 0.0) lit_monitor_set_onevalue(12, v12);
if(v13 != 0.0) lit_monitor_set_onevalue(13, v13);
if(v14 != 0.0) lit_monitor_set_onevalue(14, v14);
if(v15 != 0.0) lit_monitor_set_onevalue(15, v15);
if(v16 != 0.0) lit_monitor_set_onevalue(16, v16);
if(v17 != 0.0) lit_monitor_set_onevalue(17, v17);
if(v18 != 0.0) lit_monitor_set_onevalue(18, v18);
if(v19 != 0.0) lit_monitor_set_onevalue(19, v19);*/
}
void monitor::lit_monitor_pre_init(){
std::ofstream myfile;
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<litparam>ppar(new litparam);
lit_nfiles = 0;
gfile = "monitor/lit_log"; // check if this is the current directory
if(pp->ierr==0) lit_mfiles = new lit_mfile_t[lit_nfiles_max];
if(pp->ierr != 0) pp->litError("lit_monitor_pre_init","allocation error for lit_mfiles");
if(pp->myrank==0){
	boost::filesystem::create_directories("monitor");
	lit_file_log = pgv->lit_iopen();
	myfile.open(ppar->trim(gfile)); // check format
}
myfile.close();
}

void monitor::lit_monitor_post_init(){
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<litparam>ppar(new litparam);
std::string col(col_len, ' ');
std::string header(4096, ' ');
std::string filename(64,' '), buffer(64, ' ');
int offset;
std::fstream myfile;
bool twoLines;
if(pp->myrank != 0) return;
for(auto i=0; i<lit_nfiles; i++){
	twoLines = false;
	gfile = "monitor/"+lit_mfiles[i].filename;
	std::cout<<	gfile <<std::endl;
	lit_mfiles[i].iunit = pgv->lit_iopen();
	myfile.open(gfile);
	switch(lit_mfiles[i].freq){
		case 1:
			header.insert(0, "Step"); // format
			header.insert(1*col_len, "Time"); // format
			offset = 2*col_len;
			break;
		case 2:
			header.insert(0, "Step"); // format
			header.insert(1*col_len, "Time"); // format
			header.insert(2*col_len, "Niter"); // format
			offset = 3*col_len;
			break;
	}
	for(auto icol=0; icol < lit_mfiles[i].ncols; icol++){
		buffer.insert(0, lit_mfiles[i].header[icol]);
		auto index1 = buffer.find(' ');
		if(index1 > 0 && index1<col_len-2){ //  checked
			twoLines = true;
			col.insert(0,buffer.substr(0, index1)); // format
		}
		else{
			col.insert(0,buffer); // format
		}
		header.insert(offset+icol*col_len, ppar->trim(col));
	}
	myfile<<ppar->trim(header);
	if(twoLines){
		header = std::string(4096, ' ');
		for(auto icol=0; icol<lit_mfiles[i].ncols; icol++){
			buffer = lit_mfiles[i].header[icol];
			auto index1 = buffer.find(' ');
			if(index1 > 0 && index1 < col_len-1){
				col = buffer.substr(index1);
			}
			else{
				col = buffer;
			}
			header.insert(offset+icol*col_len, ppar->trim(col)); // format
		}
		myfile<<ppar->trim(header); // format
	}
}
}

void monitor::lit_monitor_create_gnuplot(const int iy1,const int iy2,const int iy3,\
													  const int iy4,const int iy5,const int iy6,\
													  const int iy7,const int iy8,const int iy9,\
													  const bool scnd){
/*std::vector<int> iy(9);
std::string line(4096, ' '), c2(2, ' ');
auto n=-1;
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<litparam>ppar(new litparam);
if(iy1 != 0){
	n += 1;
	iy[n] = iy1;
}
if(iy2 != 0){
	n += 1;
	iy[n] = iy2;
}
if(iy3 != 0){
	n += 1;
	iy[n] = iy3;
}
if(iy4 != 0){
	n += 1;
	iy[n] = iy4;
}
if(iy5 != 0){
	n += 1;
	iy[n] = iy5;
}
if(iy6 != 0){
	n += 1;
	iy[n] = iy6;
}
if(iy7 != 0){
	n += 1;
	iy[n] = iy7;
}
if(iy8 != 0){
	n += 1;
	iy[n] = iy8;
}
if(iy9 != 0){
	n += 1;
	iy[n] = iy9;
}
lit_ifile = lit_nfiles - 1;
auto iunit = pgv->lit_iopen();
std::ofstream myfile;
std::string strfile = "monitor/gnuplot"+ppar->trim(lit_mfiles[lit_ifile].filename);
myfile.open(strfile);
line = "p";
for(auto i=0; i<n; i++){
	if(lit_mfiles[lit_ifile].freq == 1){
		c2 = std::to_string(iy[i] + 2);
	}
	else{
		c2 = std::to_string(iy[i] + 3);
	}
	if(i==n-1){
		line = ppar->trim(line) + ' ' + char(34)+ppar->trim(lit_mfiles[lit_ifile].filename)+char(34);
		line+= " u 2:"+c2+" w l title "+char(34) + ppar->trim(lit_mfiles[lit_ifile].header[iy[i]-1])+char(34);
	}
	else{
		line = ppar->trim(line)+ ' '+char(34)+ppar->trim(lit_mfiles[lit_ifile].filename)+char(34);
		line+= " u 2:"+c2+" w l title " + char(34)+ppar->trim(lit_mfiles[lit_ifile].header[iy[i]-1])+char(34)+",";
	}
}
myfile<<ppar->trim(line)<<std::endl;
myfile<<"pause 10"<<std::endl;
myfile<<"reread"<<std::endl;
myfile.close();
auto i = pgv->lit_iclose(iunit);
if(scnd == true){
	iunit = pgv->lit_iopen();
	myfile.open(strfile);
	line = "p";
	for(auto i=0; i<n; i++){
		if(lit_mfiles[lit_ifile-1].freq == 1){
			c2 = std::to_string(iy[i]+1);
		}
		else{
			c2 = std::to_string(iy[i]+2);
		}
		if(i == n-1){
			line = ppar->trim(line)+' '+char(34)+ppar->trim(lit_mfiles[lit_ifile].filename)+char(34);
			line+=" u2 :"+c2+" w l title "+char(34)+ppar->trim(lit_mfiles[lit_ifile].header[iy[i]-2])+char(34);
		}
		else{
			line=ppar->trim(line)+' '+char(34)+ppar->trim(lit_mfiles[lit_ifile].filename)+char(34);
			line+=" u 2:"+c2+ " w l title "+char(34)+ppar->trim(lit_mfiles[lit_ifile].header[iy[i]-2])+char(34)+",";
		}
	}
	myfile<<ppar->trim(line)<<std::endl;
	myfile<<"pause 10"<<std::endl;
	myfile<<"reread"<<std::endl;
	myfile.close();
	i = pgv->lit_iclose(iunit);
}*/
}

void monitor::lit_monitor_create_file_step(const std::string filename, const int ncols){
/*std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<global_variable>pgv(new global_variable);
std::unique_ptr<litparam>ppar(new litparam);
lit_nfiles += 1;
if(lit_nfiles > lit_nfiles_max){
	pp->litError("monitor_create_file_step","too many files to monitor");
}
lit_ifile = lit_nfiles;
lit_mfiles[lit_ifile-1].filename = filename;
lit_mfiles[lit_ifile-1].freq = 1;
lit_mfiles[lit_ifile-1].ncols = ncols;
if(pp->ierr==0){
	lit_mfiles[lit_ifile-1].header = new std::string[ncols];
	lit_mfiles[lit_ifile-1].col_type = new char[ncols];
	lit_mfiles[lit_ifile-1].val = new double[ncols];
}
else	pp->litError("lit_monitor_create_file_step","allocation error for lit_mfiles");*/
}
void monitor::lit_monitor_create_file_iter(const std::string filename, const int ncols,\
														 const int ntype){
/*std::unique_ptr<parallel>pp(new parallel);
if(ntype <= 1){
	pp->litError("lit_monitor_create_file_iter","iter monitor file must be ntype>1");
}
lit_nfiles += 1;
if(lit_nfiles > lit_nfiles_max){
	pp->litError("lit_monitor_create_file_iter","too many files to monitor");
}
lit_ifile = lit_nfiles;
lit_mfiles[lit_ifile-1].filename = filename;
lit_mfiles[lit_ifile-1].freq = ntype;
lit_mfiles[lit_ifile-1].ncols= ncols;
if(pp->ierr == 0){
	lit_mfiles[lit_ifile-1].header = new std::string[ncols];
//	lit_mfiles[lit_ifile-1].col_type = std::string()
	lit_mfiles[lit_ifile-1].val = new double[ncols];
}
else
	pp->litError("lit_monitor_create_file_iter","allocation error for lit_mfiles");*/
}
void monitor::lit_monitor_select_file(std::string filename){
/*std::unique_ptr<litparam>ppar(new litparam);
std::unique_ptr<parallel>pp(new parallel);
for(lit_ifile=1; lit_ifile<lit_nfiles+1; lit_ifile++){
	if(ppar->trim(lit_mfiles[lit_ifile-1].filename) == ppar->trim(filename)) break;
}
if(lit_ifile > lit_nfiles){
	std::cout<<"filename : "<<filename<<std::endl;
	pp->litError("monitor_select_file","inknown file to monitor");
}*/
}

void monitor::lit_monitor_set_header(const int icol,const std::string header,const char col_type){
/*std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<litparam>ppar(new litparam);
if(icol > lit_mfiles[lit_ifile-1].ncols){
	std::cout<<"filename : "<<lit_mfiles[lit_ifile-1].filename<<std::endl;
	std::cout<<"column index"<<icol<<std::endl;
	pp->litError("monitor_set_header","column index too large");
}
lit_mfiles[lit_ifile-1].header[icol-1] = ppar->trim(header);
lit_mfiles[lit_ifile-1].col_type[icol-1] = col_type;*/
}

void monitor::lit_monitor_set_array_values(double *val){
//	lit_mfiles[lit_ifile-1].val = val;
}

void monitor::lit_monitor_log(const std::string text){
std::vector<int> date(8);
std::unique_ptr<parallel>pp(new parallel);
std::ofstream myfile;
myfile.open("monitor/lit_log", std::ios::app); // check if it is append mode
if(pp->myrank == 0){
	time_t rawtime;
	struct tm * timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	myfile<<asctime(timeinfo)<<std::endl;
}
}

void monitor::lit_monitor_finalize(){
std::unique_ptr<parallel>pp(new parallel);
if(pp->myrank == 0){
	for(lit_ifile=1; lit_ifile<lit_nfiles+1; lit_ifile++){
		//close a file  //use lit_ifile-1 as index do not need this part
	}
	//close a file ; do not need this part
}
}

std::string monitor::gfile = "";
double monitor::monitor_gradG_max = 0.0;
double monitor::monitor_gradG_min = 0.0;
int monitor::monitor_reinit_steps = 0;
int monitor::monitor_nSubCycles = 0;
int monitor::lit_nfiles = 0;
int monitor::lit_ifile = 0;
lit_mfile_t* monitor::lit_mfiles = nullptr;
