//Written by Salar Safarkhani

#pragma once
#include<string>
using std::string;

class litstring{
	public:
	bool string_compare(string str1, string str2);
	bool string_starts_with(string str1, string str2);
	bool compare_tokens_to_string(string tokens, string str);
	bool compare_tokens_to_tokens_old(string tokens1, string tokens2,\
												 string delimiters);
	bool compare_tokens_to_tokens(string tokens1, string tokens2);
	void convert_to_uppercase(string &str);
	private:
	inline string trim(std::string s, std::string ws=" \t\n\r") {
		if (s.empty()) return s;
		size_t a=s.find_first_not_of(ws);
		size_t b = s.find_last_not_of(ws)+1;
		return s.substr(a, b-a);
	}
   inline int counter(string str){
      string tmp = trim(str);
      return tmp.size();
   };

};

