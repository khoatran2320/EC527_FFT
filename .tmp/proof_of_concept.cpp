
#include <math.h>
#include <iostream>
#include <complex>
#include <immintrin.h>
#include <bitset>

#define PI 3.14159265358979323846
#define PERIODS 1
#define FREQ	100
using namespace std;


void proof_of_concept(){	
	int freq = 10;
	int periods = 1;
	int len = (int) 2 * PI * periods * freq;
	int index = 0;
	float input[len];
	__m256 blocks[len];

	for (int i = 0; i < len; i++) {input[i] = i;}

	for (int j = 0; j < len; j += 8) {
	//	cout << "j : " << j << endl;


	float * a_float = new (std::align_val_t(32)) float[8];
	float * b_float = new (std::align_val_t(32)) float[8];

	for (int i = 0; i < 8; i++) a_float[i] = input[index + i];//= 10;// * sin((float)i/freq);
	for (int i = 0; i < 8; i++) {
		switch(i) {
			float temp;
			case 1: 
				temp = a_float[1];
				b_float[1] = a_float[4];
				b_float[4] = temp;
				break;
			case 3: 
				temp = a_float[3];
				b_float[3] = a_float[6];
				b_float[6] = temp;
				break;
			default: b_float[i] = a_float[i];
		} 

		//cout << i << " : " << b_float[i] << endl;
		}
//	cout << i << " : " << (bs) << " : " << (bs >> 1) << endl;}// b_float[i] = 2 * input[index + i]; //sin(((float)i + index)/freq);
//	float b_float[] = {a_float[0],a_float[4], a_float[2], a_float[6], a_float[1], a_float[5], a_float[2], a_float[7]}; 

//	for (int i = 0; i < len; i++) cout << a_float[i] << endl;


	for (int size = 2; size <= 8; size *= 2) {
		size_t halfsize = size/2;
		size_t tablestep = n/size;

		for (int p = 0; p < n; p += size) {
			for (int q = p, r = 0; j < j + halfsize; j++, k += tablestep) {
				size_t l = q + halfsize; 
			}
		}
	}
	__m256 a = _mm256_load_ps(a_float);
	__m256 b = _mm256_load_ps(b_float);
	__m256 c = _mm256_sub_ps(b,a);
	__m256 d = _mm256_mul_ps(b,a);
	__m256 e = _mm256_add_ps(b,a);
	blocks[0] = c;
	float * c_float = new (std::align_val_t(32)) float[8];	


	_mm256_store_ps(c_float,c);

//	for (int i = 0; i < 8; i++) cout << c_float[i] << endl;
	index += 8;
}
	float * d_float = new (std::align_val_t(32)) float[8];	

	_mm256_store_ps(d_float, blocks[0]);
//	for (int i = 0; i < 8; i++) cout << d_float[i] << endl;
}

int main(){
	proof_of_concept();
	return 0;
}