//Written by Salar Safarkhani

#pragma once
#include<string>
class InputFile;

class Simulation{
public:
	Simulation();
	~Simulation();
	void read_input(const std::string &file);
	void command(const std::string &line);
private:
	InputFile *_input;
};
