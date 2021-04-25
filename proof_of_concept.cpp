
#include <math.h>
#include <iostream>
#include <complex>
#include <immintrin.h>

#define PI 3.14159265358979323846
#define PERIODS 1
#define FREQ	100
using namespace std;

void proof_of_concept(){	
	int freq = 10;
	int periods = 1;
	int len = (int) 2 * PI * periods * freq;
	int index = 0;

	__m256 blocks[len];

	for (int j = 0; j < len; j += 8) {
	//	cout << "j : " << j << endl;

	float * a_float = new (std::align_val_t(32)) float[8];
	float * b_float = new (std::align_val_t(32)) float[8];

	for (int i = 0; i < 8; i++) a_float[i] = 10;// * sin((float)i/freq);
	for (int i = 0; i < 8; i++) b_float[i] = sin(((float)i + index)/freq);


//	for (int i = 0; i < len; i++) cout << a_float[i] << endl;


	__m256 a = _mm256_load_ps(a_float);
	__m256 b = _mm256_load_ps(b_float);
	__m256 c = _mm256_sub_ps(b,a);

	blocks[0] = c;
	float * c_float = new (std::align_val_t(32)) float[8];	


	_mm256_store_ps(c_float,c);

	//for (int i = 0; i < 8; i++) cout << c_float[i] << endl;
	index += 8;
}
	float * d_float = new (std::align_val_t(32)) float[8];	

	_mm256_store_ps(d_float, blocks[0]);
	for (int i = 0; i < 8; i++) cout << d_float[i] << endl;
}