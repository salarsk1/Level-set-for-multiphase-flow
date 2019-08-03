//Written by Salar Safarkhani

#include"param.h"
#include<mpi.h>
#include"parallel.h"
#include<iostream>
#include<fstream>
#include<algorithm>
#include<sstream>
void litparam::init_param(int argc, char *argv[], const string& input_filename){
std::unique_ptr<parallel>pp(new parallel);
string argv0, filename;
int backslash, ierr, my_flag, flag;
first_param_ptr = nullptr;
parse_commandline_param(argc, argv, my_flag);
MPI_Allreduce(&my_flag, &flag, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD);
if(flag != 0){
	if(pp->myrank==0)
	 std::cout<<"Warning: problem parsing commandline-skipping"<<std::endl;
	destroy_param();
	first_param_ptr = nullptr;
};
if(input_filename != "\0")	filename = input_filename;
else{
	argv0 = argv[1];
	trim(argv0);
	backslash = argv0.rfind("/");
	if(backslash<0 || backslash>argv0.size()){
		filename =  trim(argv0)+".in";
	}
	else{
		std::string stmp;
		for(auto i=0; i<argv0.size(); i++){
			stmp = argv0[backslash+i];
			filename+=trim(stmp);
		}
		filename += ".in";
	}
}
read_param_file(filename);
};

void litparam::destroy_param_recursive(param_t *param){
std::unique_ptr<parallel>pp(new parallel);
if(param->next != nullptr){
	destroy_param_recursive(param->next);
	if(pp->ierr==0)
		delete param->next;
	else
		pp->litError("destroy_param_recursive","deallocation error for\
						  param->next");
};
if(pp->ierr==0)
	delete []param->argv; // checked
else
	pp->litError("destroy_param_recursive","deallocation error for\
					  param->argv");
};

void litparam::destroy_param(){
std::unique_ptr<parallel>pp(new parallel);
if(first_param_ptr != nullptr){
	destroy_param_recursive(first_param_ptr);
	if(pp->ierr==0)
		delete first_param_ptr;
	else
	 pp->litError("destroy_param","deallocation error for first_param_ptr");
};
};

void litparam::parse_commandline_param(int argc, char * argv[], int &ierr){
int iarg, len_trim_token, iargc;
int iargv[MAX_TOKEN_COUNT+1];
string token, line;
#ifdef NO_GETARG
	std::cout<<"Warning: NO_GETARG is defined: skipping parse_commandline_param."<<std::endl;
	return;
#else
	ierr  = 0;
	iarg  = 1;
	token = argv[iarg];
	trim(token);
	len_trim_token = token.size();
	iargc = 0;
	iargv[0] = 1;
	while(len_trim_token > 0){
		if(token[0] == '-' && token[1] == '-'){
			if(len_trim_token < 3){
				ierr = -1;
				return;
			}
			if(iargc > 0){
				trim(line); // checked
				if((line[0]=='i'&&line[1]=='n'&&line[2]=='c'&&line[3]=='l'&&line[4]=='u'&&\
					line[5]=='d'&&line[6]=='e')||\
					line[0]=='I'&&line[1]=='N'&&line[2]=='C'&&line[3]=='L'&&line[4]=='U'&&\
					line[5]=='D'&&line[6]=='E'){
					if(iargc != 2){
						ierr = -1;
						return;
					}
					int siz = iargv[2]-2-iargv[1]+1;
					std::string temp_line = line.substr(iargv[1]-1,siz);
					read_param_file(temp_line);
				}
				else{
					add_param_line(line, iargc, iargv);
				}
				iargc = 0;
			}
			std::string token_temp = &token[2];
			trim(token_temp);
			line = token_temp;
			iargc += 1;
			iargv[iargc] = line.size()+2;
		}
		else{
			if(iargc == 0){
				ierr = -1;
				return;
			}
			trim(line);
			trim(token);
			line = line + " " + token;
		}
		iarg += 1;
		token = argv[iarg];
		std::string token_temp = token;
		trim(token_temp);
		len_trim_token = token_temp.size();
	}
	if(iargc > 0){
		trim(line); // checked
		if((line[0]=='i'&&line[1]=='n'&&line[2]=='c'&&line[3]=='l'&&line[4]=='u'&&\
			line[5]=='d'&&line[6]=='e')||\
			line[0]=='I'&&line[1]=='N'&&line[2]=='C'&&line[3]=='L'&&line[4]=='U'&&\
			line[5]=='D'&&line[6]=='E'){
			if(iargc != 2){
				ierr = -1;
				return;
			}
			std::string line_temp;
			line_temp.assign(line, iargv[1]-1, iargv[2]-iargv[1]-1);
			read_param_file(line_temp);
		}
		else{
			add_param_line(line, iargc, iargv);
		}
	}
#endif
};

param_t* litparam::get_new_param(){
std::unique_ptr<parallel>pp(new parallel);
param_t *param;
if(first_param_ptr == nullptr){
	if(pp->ierr==0) first_param_ptr = new param_t;
	else pp->litError("get_new_param", "allocation error for first_param_ptr");
	param = first_param_ptr;
}
else{
	param = first_param_ptr;
	while(param->next != nullptr) param = param->next;
	if(pp->ierr==0) param->next = new param_t;
	else pp->litError("get_new_param","allocation error for param->next");
	param = param->next;
}
param->next = nullptr;
param->line = "";
param->argc = 0;
param->request_count = 0;
param->argv = nullptr;
return param;
};

void litparam::build_mpi_filename(std::string &mpi_filename, std::string filename){
#ifdef MPI_UFS_PREFIX
	trim(filename);
	mpi_filename = "ufs:"+filename;
#else
	trim(filename);
	mpi_filename = filename;
#endif
};

int litparam::get_next_token(int pos, int max_pos, int& IACHAR_TAB, \
						 int& IACHAR_RETURN, int& IACHAR_SPACE, int& IACHAR_COMMA, \
						 int& IACHAR_EQUALS, int& IACHAR_SEMICOLON, int& IACHAR_HASH, \
						 int& IACHAR_QUOTE, string& token, bool& comment_mode, string char_buffer,\
						 MPI_File fh, MPI_Status status, int ierr){
std::unique_ptr<parallel> pp(new parallel);
int flag = 0;
int token_pos = 0;
bool quote_mode = false;
int ic;
std::fill(token.begin(), token.end(), ' ');
while(true){
	while(pos < max_pos){
		ic = (int)(char_buffer[pos]);
		if(quote_mode){
			pos += 1;
			token_pos += 1;
			if(token_pos >= MAX_TOKEN_LEN){
				std::cout<<"Error: increase MAX_TOKEN_LEN"<<std::endl;
				pp->parallel_kill(0);
			}
			token[token_pos-1] = char_buffer[pos-1]; // checked
			if(ic == IACHAR_QUOTE) return flag;
		}
		else if(comment_mode){
			pos += 1;
			if(ic == IACHAR_RETURN) comment_mode = false;
		}
		else{
			if(ic == IACHAR_QUOTE){
				if(token_pos != 0){
					std::cout<<"Error: quote not at token_pos == 0"<<std::endl;
					pp->parallel_kill(0);
				}
				pos += 1;
				token_pos = 1;
				token[token_pos-1] = '"';
				quote_mode = true;
			}
			else if(ic == IACHAR_SPACE||ic==IACHAR_COMMA||ic==IACHAR_TAB||ic==IACHAR_EQUALS){
				if(token_pos > 0){
					pos += 1;
					return flag;
				}
				else{
					pos += 1;
				}
			}
			else if(ic==IACHAR_RETURN||ic==IACHAR_SEMICOLON){
				if(token_pos > 0){
					pos += 1;
					flag=  1;
					return flag;
				}
				else{
					pos += 1;
					flag = 2;
					return flag;
				}
			}
			else if(ic==IACHAR_HASH){
				comment_mode = true;
				if(token_pos > 0){
					pos += 1;
					flag = 1;
					return flag;
				}
				else{
					pos += 1;
					flag = 2;
					return flag;
				}
			}
			else{
				pos += 1;
				token_pos += 1;
				if(token_pos >= MAX_TOKEN_LEN){
					std::cout<<"Error: increase MAX_TOKEN_LEN."<<std::endl;
					pp->parallel_kill(0);
				}
				token[token_pos-1] = char_buffer[pos-1];
			}
		}
	}
	ierr = MPI_File_read(fh, &char_buffer[0], char_buffer.size(), MPI_CHARACTER, &status);
	ierr = MPI_Get_count(&status, MPI_CHARACTER, &max_pos);
	pos = 0;
	if(max_pos == 0){
		if(token_pos > 0){
			flag = 1;
			return flag;
		}
		else{
			flag = -1;
			return flag;
		}
	}
}
};

void litparam::read_param_file(string filename){
size_t CHAR_BUFFER_LEN = 1024;
int IACHAR_TAB = 9;
int IACHAR_RETURN = 10;
int IACHAR_SPACE = (int) ' ';
int IACHAR_COMMA = (int) ',';
int IACHAR_EQUALS = (int) '=';
int IACHAR_SEMICOLON = (int) ';';
int IACHAR_HASH = (int) '#';
int IACHAR_QUOTE = (int) '"';
std::unique_ptr<parallel>pp(new parallel);
MPI_File fh;
int ierr, pos, max_pos;
MPI_Status status; // checked
string token;
string char_buffer;
string line;
int argc;
int argv[MAX_TOKEN_COUNT+1];
string mpi_filename;
bool comment_mode;
if(pp->myrank == 0) std::cout<<"in_read_param_file "<<trim(filename)<<"..."<<std::endl;
recursion_level = recursion_level + 1;
if(recursion_level > MAX_RECURSION_LEVEL){
	std::cout<<"Error: increase MAX_RECURSION_LEVEL, or check includes."<<std::endl;
	pp->parallel_kill(0);
}
if(filename[0] == '"'){
	trim(filename);
	ierr = trim(filename).size();
	if(filename[ierr-1] != '"'){
		std::cout<<"Error: expected double quote at end: "<<filename<<std::endl;
		pp->parallel_kill(0);
	}
	build_mpi_filename(mpi_filename, &filename[1]);
}
else{
	build_mpi_filename(mpi_filename, filename);
}
 ierr=MPI_File_open(MPI_COMM_WORLD,&mpi_filename[0],MPI_MODE_RDONLY,pp->MPI_IO_LIT_FILE_HINT,&fh);
if(ierr!=0){
	if(pp->myrank==0) std::cout<<"Could not find/open file. skipping."<<std::endl;
	return;
}
pos = 0;
max_pos = 0;
argc = 0;
argv[0] = 1;
comment_mode = false;
while(true){
	int sm = get_next_token(pos,max_pos,IACHAR_TAB,IACHAR_RETURN,\
									IACHAR_SPACE,IACHAR_COMMA,IACHAR_EQUALS,IACHAR_SEMICOLON,\
									IACHAR_HASH,IACHAR_QUOTE,token,comment_mode,char_buffer,\
                 				fh,status,ierr);
	if(sm == -1){
		break;
	}
	else if(sm == 0){
		if(argc == 0){
			trim(token);
			line = token;
		}
		else{
			trim(token);
			trim(line);
			line = line+" "+token;
		}
		argc += 1;
		argv[argc] = line.size()+2;
	}
	else if(sm == 1){
		if(argc == 0){
			trim(token);
			line = token;
		}
		else{
			trim(token);
			trim(line);
			line = line + " " + token;
		}
		argc += 1;
		argv[argc] = line.size() + 2;
      if((line[0]=='i'&&line[1]=='n'&&line[2]=='c'&&line[3]=='l'&&line[4]=='u'&&\
	       line[5]=='d'&&line[6]=='e')||\
     		 line[0]=='I'&&line[1]=='N'&&line[2]=='C'&&line[3]=='L'&&line[4]=='U'&&\
	       line[5]=='D'&&line[6]=='E'){
			if(argc != 2){
				std::cout<<"Error: expect only one filename with include."<<std::endl;
				pp->parallel_kill(0);
			}
			string line_temp;
			line_temp.assign(line, argv[1]-1, argv[2]-argv[1]-1);
		}
		else add_param_line(line, argc, argv);
		argc = 0;
	}
	else if(sm == 2){
		if(argc > 0){
			if((line[0]=='i'&&line[1]=='n'&&line[2]=='c'&&line[3]=='l'&&line[4]=='u'&&\
	         line[5]=='d'&&line[6]=='e')||\
   	      line[0]=='I'&&line[1]=='N'&&line[2]=='C'&&line[3]=='L'&&line[4]=='U'&&\
      	   line[5]=='D'&&line[6]=='E'){
				if(argc != 2){
					std::cout<<"Error: expect only one filename with include."<<std::endl;
					pp->parallel_kill(0);
				}
           	string line_temp; 
           	line_temp.assign(line, argv[1]-1, argv[2]-argv[1]-1);
			}
			else add_param_line(line, argc, argv);
			argc = 0;
		}
	}
}
ierr = MPI_File_close(&fh);
recursion_level -= 1;
};
void litparam::add_param_line(const string line,int argc,int* argv){
std::unique_ptr<parallel>pp(new parallel);
param_t* param;
param = get_new_param();
param->line = line;
param->argc = argc;
if(pp->ierr==0) param->argv = new int[argc+1];
else pp->litError("add_param_line","allocation error for param->argv");
for(auto i=0; i<argc+1; i++) param->argv[i] = argv[i];
};

param_t* litparam::get_param(const string name){
param_t* param;
param = first_param_ptr;
string s1 = trim(name);
while(param != nullptr){
	if(param->line.find(name) == 0 && s1.size()==param->argv[1]-2){
		param->request_count += 1;
		return param;
	}
	param = param->next;
}
return param;
};

void litparam::set_next_param(param_t* param, const string name){
if(param == nullptr) return;
param = param->next;
string s1 = trim(name);
while(param != nullptr){
	if(param->line.find(name) == 0 && s1.size() == param->argv[1]-2){
		param->request_count += 1;
		break;
	}
	param = param->next;
}
};
int litparam::get_param_type(const param_t &param, const int iarg){
int param_type, ierr;
int argv_f,argv_l;
double test_real_wp;
std::ifstream myfile;
if(iarg >= param.argc){
	param_type = NO_PARAM_TYPE;
	return param_type;
}
argv_f = param.argv[iarg];
argv_l = param.argv[iarg+1] - 2;
string st = param.line.substr(argv_f-1);
myfile.open(st); // checked
if(myfile.good()){
	myfile>>test_real_wp; // checked
}
if(myfile.good() && int(param.line[argv_f-1])<65){
	param_type = NUMERICAL_PARAM_TYPE;
	return param_type;
}
param_type = STRING_PARAM_TYPE;
return param_type;
};

void litparam::dump_param(){
param_t* param;
std::cout<<"*************** params ****************"<<std::endl;
param = first_param_ptr;
while(param != nullptr){
	std::cout<<trim(param->line)<<std::endl;
	param = param->next;
}
};

void litparam::dump_param_usage(){
int ierr, pcount;
std::unique_ptr<parallel> pp(new parallel);
param_t* param;
std::vector<int> my_used, used;
pcount = 0;
param = first_param_ptr;
while(param != nullptr){
	pcount += 1;
	param = param->next;
}
if(pcount > 0){
	if(ierr == 0) my_used.resize(pcount);
	else if(ierr != 0) pp->litError("dump_param_usage","allcoation error for my_used");
	for(auto i=0; i<pcount; i++) my_used[i] = 0;
	pcount = 0;
	param = first_param_ptr;
	while(param != nullptr){
		pcount += 1;
		my_used[pcount - 1] = param->request_count;
		param = param->next;
	}
	if(pp->myrank == 0){
		if(ierr==0)	used.resize(pcount);
		else pp->litError("dump_param_usage", "allocation error for used");
	}
	ierr = MPI_Reduce(&my_used[0], &used[0], pcount, MPI_INTEGER, MPI_MAX,0, MPI_COMM_WORLD);
}
if(pp->myrank==0){
	std::cout<<"*********** param usage *************"<<std::endl;
	pcount = 0;
	param = first_param_ptr;
	while(param != nullptr){
		pcount += 1;
		int siz = param->argv[1]-param->argv[0]-1;
		string st = param->line.substr(param->argv[0]-1,siz);
		std::cout<<trim(st)<<", max request_count= "<<used[pcount-1]<<std::endl;
		param = param->next;
	}
	std::cout<<"************ end of param usage *************"<<std::endl;
	if(pcount > 0){
		if(ierr == 0){
			used.clear();
			used.shrink_to_fit();
		}
		else pp->litError("dump_param_usage","deallocation error for used");
	}
}
if(pcount > 0){
	if(ierr == 0){
		my_used.clear();
		my_used.shrink_to_fit();
	}
	else pp->litError("dump_param_usage","deallocation error for my_used");
}
};

bool litparam::check_param(const string name){
std::unique_ptr<litstring>pstr(new litstring);
param_t* param;
param = first_param_ptr;
while(param != nullptr){
	string st = param->line.substr(0, param->argv[1]-2);
	if(pstr->string_compare(st, name)){
		param->request_count += 1;
		return true;
	}
	param = param->next;
}
return false;
};

bool litparam::get_logical_param(const string name, bool defa, int arg){
std::unique_ptr<litstring>pstr(new litstring);
std::unique_ptr<parallel> pp(new parallel);
param_t* param;
param = first_param_ptr;
bool ret;
int ierr;
std::ifstream myfile;
std::string st2, st3;
while(param != nullptr){
	string st1 = param->line.substr(0,param->argv[1]-2);
	if(pstr->string_compare(st1, name)){
		if(param->argc==1) ret = true;
		else{
			st2 = param->line.substr(param->argv[1]-1);
			if(pp->ierr==0){
				st3 = lower(st2);
				std::istringstream(st3)>>std::boolalpha>>ret;
			}
			else{
				std::cout<<"Format for param "<<trim(name)<<"does not match requested logical";
				std::cout<<std::endl;
				pp->parallel_kill(0);
			}
		}
		param->request_count += 1;
		return ret;
	}
	param = param->next;
}
if(arg != 1) return defa;
else return false;
};

int litparam::get_integer_param(const string name, int defa, int arg){
int ierr, ret;
std::unique_ptr<parallel>pp(new parallel);
std::unique_ptr<litstring>pstr(new litstring);
param_t* param;
param = first_param_ptr;
std::string st3;
while(param != nullptr){
	string st = param->line.substr(0, param->argv[1]-2);
	if(pstr->string_compare(st, name)){
		if(param->argc == 1){
			std::cout<<"Error: param "<<trim(name)<<" has no associated value"<<std::endl;
			pp->parallel_kill(0);
		}
		st3 = param->line.substr(param->argv[1]-1);
		if(pp->ierr==0) ret = std::stoi(st3);
		else{
			std::cout<<"Format for param "<<trim(name)<<"does not match requested integer"<<std::endl;
			pp->parallel_kill(0);
		}
		param->request_count += 1;
		return ret; // checked
	}
	param = param->next;
}
if(arg > 1) return defa;
else{
	std::cout<<"Error: could not find param "<<name<<std::endl;
	pp->parallel_kill(0);
	return 1;
}
};

double litparam::get_real_param(const string name, double defa, int arg){
int ierr;
std::unique_ptr<litstring>pstr(new litstring);
std::unique_ptr<parallel>pp(new parallel);
param_t* param;
param = first_param_ptr;
double ret;
std::string st4;
while(param != nullptr){
	string st = param->line.substr(0, param->argv[1]-3);
	if(pstr->string_compare(st, name)){
		if(param->argc == 1){
			std::cout<<"Error: param "<<trim(name)<<" does not associated value"<<std::endl;
			pp->parallel_kill(0);
		}
		st4 = param->line.substr(param->argv[1]-1);
		if(pp->ierr==0)ret = std::stod(st4);
		else{
			std::cout<<"Format for param "<<trim(name)<<" does not match requested real."<<std::endl;
			pp->parallel_kill(0);
		}
		param->request_count += 1;
		return ret;
	}
	param = param->next;
}
if(arg >1) return defa;
else{
	std::cout<<"Error: could not find real param "<<name<<std::endl;
	pp->parallel_kill(0);
}
};

string litparam::get_string_param(string name, string defa, int arg){
string str;
int ierr;
std::unique_ptr<litstring>pstr(new litstring);
std::unique_ptr<parallel>pp(new parallel);
param_t* param;
param=first_param_ptr;
std::ifstream myfile;
while(param != nullptr){
	string st = param->line.substr(0, param->argv[1]-2);
	if(pstr->string_compare(st, name)){
		std::cout<<"Error: param "<<trim(name)<<" has no associated value"<<std::endl;
		pp->parallel_kill(0);
		string st1 = param->line.substr(param->argv[1]-1);
		myfile.open(st1);
		myfile>>str;
		if(!myfile.good()){
			std::cout<<"Format is not appropriate"<<std::endl;
			pp->parallel_kill(0);
		}
		return str;
	}
	param = param->next;
}
if(arg>1) return trim(defa);
else{
	std::cout<<"Error: could not find string param "<<name<<std::endl;
	pp->parallel_kill(0);
}
};

std::array<int, 3> litparam::get_integer3_param(string name, std::vector<int> defa, int arg){
std::array<int, 3> iv3;
int ierr;
std::unique_ptr<litstring>pstr(new litstring);
std::unique_ptr<parallel>pp(new parallel);
param_t* param;
std::ifstream myfile;
param = first_param_ptr;
while(param != nullptr){
	string st1 = param->line.substr(0, param->argv[1]-2);
	if(pstr->string_compare(st1, name)){
		if(param->argc == 1){
			std::cout<<"Error: param "<<trim(name)<<" has no associated value."<<std::endl;
			pp->parallel_kill(0);
		}
		string st2 = param->line.substr(param->argv[1]-1);
		myfile.open(st2);
		for(auto &i:iv3) myfile>>i;
		if(!myfile.good()){
			std::cout<<"Format for param "<<trim(name)<<" does not match requested 3D integer."<<std::endl;
			pp->parallel_kill(0);
		}
		param->request_count += 1;
		return iv3;
	}
	param = param->next;
}
if(arg > 1){
	for(auto i=0; i<3; i++) iv3[i] = defa[i];
}
else{
	std::cout<<"Error: could not find param "<<name<< std::endl;
	pp->parallel_kill(0);
}
};

std::array<double, 3> litparam::get_real3_param(const string  name,\
										const  std::vector<double> defa,int arg){
std::array<double, 3> rv3;
int ierr;
std::unique_ptr<litstring>pstr(new litstring);
std::unique_ptr<parallel>pp(new parallel);
param_t* param;
std::ifstream myfile;
param = first_param_ptr;
string st1, strtemp;
string::size_type sz, sz1;
if(first_param_ptr == nullptr) std::cout<<"we have null ptr"<<std::endl;
while(param != nullptr){
	st1 = param->line.substr(0, param->argv[1]-2);
	if(pstr->string_compare(st1, name)){
		if(param->argc == 1){
			std::cout<<"Error: param "<<trim(name)<<" has no associated value"<<std::endl;
			pp->parallel_kill(0);
		}
		strtemp = param->line.substr(param->argv[1]-1);
      if(pp->ierr==0){
         string strtemp = param->line.substr(param->argv[2]-1);
         rv3[0] = std::stod(strtemp, &sz);
         rv3[1] = std::stod(strtemp.substr(sz), &sz1);
         rv3[2] = std::stod(strtemp.substr(sz+sz1));
      }
		else{
		 std::cout<<"Format for param "<<trim(name)<<" does not match requested 3D double"<<std::endl;
		 pp->parallel_kill(0);
		}
		param->request_count += 1;
		return rv3;
	}
	param = param->next;
}
if(arg > 1){
	for(auto i=0; i<3; i++) rv3[i] = defa[i];
}
else{
	std::cout<<"Error: could not find param "<<name<< std::endl;
	pp->parallel_kill(0);
}
};

void litparam::set_param(string &value, const string name, int index,const string defa){
std::unique_ptr<litstring>pstr(new litstring);
std::unique_ptr<parallel>pp(new parallel);
param_t* param;
int index_copy;
string st1, st2;
std::ifstream myfilel;
if(index != -1) index_copy = index;
else index_copy = 2;
param = first_param_ptr;
while(param != nullptr){
	st1 = param->line.substr(0, param->argv[1]-2);
	if(pstr->string_compare(st1, name)){
		if(param->argc < index_copy){
			std::cout<<"Error: param "<<trim(name)<<" has no value associated with index: "<<index_copy<<std::endl;
			pp->parallel_kill(0);
		}
		int siz = param->argv[index_copy]-param->argv[index_copy-1]-1;
		value = param->line.substr(param->argv[index_copy-1],siz);
		param->request_count += 1;
		return;
	}
	param = param->next;
}
if(defa != "\0") value = defa;
else{
	std::cout<<"Error: could not find param "<<name<<std::endl;
	pp->parallel_kill(0);
}
};

void litparam::set_param_token(string &value, const param_t& param, int index){
int index_copy;
std::unique_ptr<parallel>pp(new parallel);
if(index != -1) index_copy = index;
else index_copy = 2;
if(param.argc < index_copy){
	std::cout<<"Error: param "<<trim(param.line)<<" has no value associated with index: "<<index_copy<<std::endl;
	pp->parallel_kill(0);
}
else{
	int siz = param.argv[index_copy]-param.argv[index_copy-1]-1;
	value = param.line.substr(param.argv[index_copy-1], siz);
	return;
}
};
param_t* litparam::first_param_ptr = nullptr;
