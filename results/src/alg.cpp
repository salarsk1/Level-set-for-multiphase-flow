#include"alg.h"
double alg::minval(int size, std::vector<std::vector<double>>& v){
//	auto s1 = v.size();
//	auto s2 = v[0].size();
//	std::vector<double> v_temp(s1*s2);
//	for(auto i=0; i<s1; i++){
//		for(auto j=0; j<s2; j++){
//			v_temp[i+j*s1] = v[i][j];
//		}
//	}
//   double res = std::numeric_limits<double>::max();
//	for(auto i=0; i<s1*s2; i++){
//      res = std::min(res, std::abs(v_temp[i]));
//	}
   std::vector<double>::const_iterator it;
   double res = std::numeric_limits<double>::max();
   for(auto i=0; i<size; i++){
      it = std::min_element(v[i].begin(), v[i].end());
      res = std::min(res, std::abs(*it));
   }
   double w = res;
   res = std::numeric_limits<double>::max();
   for(auto i=0; i<size; i++){
      it = std::max_element(v[i].begin(), v[i].end());
      res = std::min(res, std::abs(*it));
   }
   w = std::min(w, res);
	return w;
}
double alg::maxval(int size, std::vector<std::vector<double>>& v){
//	auto s1 = v.size();
//	auto s2 = v[0].size();
//	std::vector<double> v_temp(s1*s2);
//	for(auto i=0; i<s1; i++){
//		for(auto j=0; j<s2; j++){
//			v_temp[i+j*s1] = v[i][j];
//		}
//	}
//   double res = std::numeric_limits<double>::min();
//	for(auto i=0; i<s1*s2; i++){
//      res = std::max(res, std::abs(v_temp[i]));
//	}
   std::vector<double>::const_iterator it;
   double res = std::numeric_limits<double>::min();
   for(auto i=0; i<size; i++){
      it = std::max_element(v[i].begin(), v[i].end());
      res = std::max(res, std::abs(*it));
   }
   double w = res;
   res = std::numeric_limits<double>::min();
   for(auto i=0; i<size; i++){
      it = std::min_element(v[i].begin(), v[i].end());
      res = std::max(res, std::abs(*it));
   }
   w = std::max(w, res);
	return w;
}

