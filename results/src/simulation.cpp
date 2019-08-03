//Written by Salar Safarkhani

#include"simulation.h"
#include"param.h"
#include"test/test_suit.h"
#include<exception>
#include<iostream>
#include<fstream>
using namespace std;

class InputFile{
public:
	string filename;
	int current_line_number;
	bool echo;
};

Simulation::Simulation(): _input(new InputFile){
	_input->echo = false;
}

Simulation::~Simulation() {
	try {
		delete _input;
	}
	catch (string &s) {
	cout << "Error clearing memory.\n";
	cout << s << "\n";
	}
}
// Reads the simulation input file, and runs commands one line at a time.
// Anything after a # in a line is a comment.
void Simulation::read_input(const string &file) {
	std::unique_ptr<litparam> pparam(new litparam);
	_input->filename = file;
	_input->current_line_number = 0;
	fstream fid(file.c_str(), ios::in);
	if (!fid) {
		cout << "Error: cannot open input file: " << file << "\n";
	}
	while (fid) {
		// Reads next line, strips any comments, and converts to lower case.
		string line = pparam->trim(pparam->lower(pparam->read_line(fid)));
		_input->current_line_number++;
		if (line.empty()) continue;
		if (_input->echo) cout << "READ COMMAND: " << line << "\n";
		// Line is not empty - run it.
		try {
			if (line == "exit" || line == "quit") break;
			command(line);
		}
		catch (const char *msg) {
		cout << "Error in command: \n";
		cout << "  " << line << "\n" << msg << "\n";
		exit(1);
		}
		catch (string &msg) {
			if (_input->filename.size()) {
				cout << "Error in line " << _input->current_line_number;
				cout << " of " << _input->filename << "\n";
			}
			else cout << "Error in command: \n";
			cout << "  " << line << "\n" << msg << "\n";
			exit(1);
		}
	}
}
// Processes a single line.
void Simulation::command(const string &line) {
	unique_ptr<litparam> pparam(new litparam);
	string command = pparam->split(line)[0];
	if (command == "test")	lit_tests::cmd_test(line);
	else throw string("Unknown command: " + line);
}

