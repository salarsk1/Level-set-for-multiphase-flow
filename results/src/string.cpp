//Written by Salar Safarkhani

#include<algorithm>
#include"string.h"
bool litstring::string_compare(string str1, string str2){
auto s1 = str1.find(trim(str2));
auto s2 = str2.find(trim(str1));
if(s1==0 && s2==0)
	return true;
else
	return false;
};
bool litstring::string_starts_with(string str1, string str2){
if(str1.find(trim(str2))==0)
	return true;
else
	return false;
};
bool litstring::compare_tokens_to_string(string tokens, string str){
int len_trim_tokens, istart, iend;
len_trim_tokens = counter(tokens);
iend = 0;
do {
	for(istart = iend; istart<len_trim_tokens; istart++){
		string s1 = ":";
		string s2;
		s2 += tokens[istart];
		if(s1.find(s2)!=0)
			break;
	};
	for(iend = istart; iend<len_trim_tokens; iend++){
		string s1 = ":";
		string s2;
		s2 += tokens[iend];
		if(s1.find(s2)==0)
			break;
	};
	if(istart < iend){
		if(string_compare(tokens.substr(istart-1, iend-istart), str)){
			return true;
		};
	}
	else{
		break;
	};
}
while(1); // checked
return false;
};

bool litstring::compare_tokens_to_tokens_old(string tokens1,string\
															tokens2, string delimiters){
int len_trim_tokens1, istart1, iend1;
int len_trim_tokens2, istart2, iend2;
len_trim_tokens1 = counter(tokens1);
len_trim_tokens2 = counter(tokens2);
iend1 = 0;
iend2 = 0;
do{
	for(istart1=iend1;istart1<len_trim_tokens1; istart1++){
		string s1;
		s1 += tokens1[istart1];
		if(delimiters.find(trim(s1))>delimiters.size() || \
			delimiters.find(trim(s1))<0) break;
	}			
	for(iend1=istart1;iend1<len_trim_tokens1; iend1++){
		string s1;
		s1 += tokens1[iend1];
		if(delimiters.find(trim(s1))<delimiters.size() && \
			delimiters.find(trim(s1))>0) break;
	}
	if(istart1<iend1){
		iend2 = 0;
		do{
			for(istart2=iend2; istart2<len_trim_tokens2; istart2++ ){
				string s1;
				s1 += tokens2[iend2];
				if(delimiters.find(trim(s1))<0 || delimiters.find(trim(s1))>\
					delimiters.size()) break;
			}
			for(iend2=istart2; iend2<len_trim_tokens2; iend2++ ){
				string s1;
				s1 += tokens2[iend2];
				if(delimiters.find(trim(s1))>0 || delimiters.find(trim(s1))<\
					delimiters.size()) break;
			}
			if(istart2<iend2){
				string s1;
				string s2;
				for(auto i=istart1; i<iend1-1;i++) s1 += tokens1[i];
				for(auto i=istart2; i<iend2-1;i++) s1 += tokens2[i];
				if(string_compare(s1,s2)) return true;
			}
			else{
				break;
			};

		}
		while(1);

	}
	else{
		break;
	};
}
while(1);
return false;
};

bool litstring::compare_tokens_to_tokens(string tokens1, string tokens2){
auto len_trim_tokens1 = counter(tokens1);
auto len_trim_tokens2 = counter(tokens2);
int  istart1, iend1, istart2, iend2;
iend1 = 0;
do{
	for(istart1=iend1; istart1<len_trim_tokens1;istart1++){
		string s1 = ":";
		string s2 = "";
		s2 += tokens1[istart1];
		if(s1.find(s2) !=0) break;
	}
	for(iend1=istart1; iend1<len_trim_tokens1;iend1++){
		string s1 = ":";
		string s2 = "";
		s2 += tokens1[iend1];
		if(s1.find(s2) ==0) break;
	}
	if(istart1<iend1){
		iend2 = 0;
		do{
			for(istart2=iend2;istart2<len_trim_tokens2;istart2++){
				string s1 = ":";
				string s2 = "";
				s2 += tokens2[istart2];
				if(s1.find(s2)!=0) break;
			}
			for(iend2=istart2;iend2<len_trim_tokens2;iend2++){
				string s1 = ":";
				string s2 = "";
				s2 += tokens2[iend2];
				if(s1.find(s2)==0) break;
			}
			if(istart2<iend2){
				string s1 = "";
				string s2 = "";
				for(auto i=istart1;i<iend1-1; i++) s1 += tokens1[i];
				for(auto i=istart2;i<iend2-1; i++) s2 += tokens2[i];
				if(string_compare(s1,s2)) return true;
			}
			else
				break;
		}
		while(1);
	}
	else
		break;
}
while(1);
return false;
};

void litstring::convert_to_uppercase(string &str){
std::transform(str.begin(), str.end(), str.begin(),\
			 static_cast<int(*)(int)>(toupper));



};
