//Written by Salar Safarkhani

#pragma once
template<class TY>
class arrayND{
	public:
		inline void getsize(TY ***typtr, size_t &d1, size_t &d2, size_t &d3){
			d1 = sizeof(typtr) / sizeof(typtr[0]);
			d2 = sizeof(typtr[0]) / sizeof(typtr[0][0]);
			d3 = sizeof(typtr[0][0]) / sizeof(TY);
		}
  
};
