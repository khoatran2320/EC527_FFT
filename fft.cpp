#include <math.h>
#include <iostream>
#include <complex>
#include <immintrin.h>

#define PI 3.14159265358979323846
#define PERIODS 1
#define FREQ	100
using namespace std;
// g++ fft.cpp -o fft_cpp -w -lm -mavx2 -std=c++17 && ./fft_cpp 
// https://www.alcf.anl.gov/files/ken_intel_compiler_optimization.pdf

//__m256d generate_sine_wave(float periods, float freq){
//	// X[samples] =
//
//	__mm128 X[(int) (2 * PI * periods * freq) + 1] ;// __attribute__((aligned(64)));	
//	__mm256 x __mm_set_ss(1.0);
//	__mm256 sine_wave
//	//__assume_aligned(X,64);
	//double samples = 2 * PI * periods * freq;
	//for (int i = 0; i < (int) samples; i++) {
	//	double x = (double) (i/freq);
	//	double y =  sin(x);
//
//	//	//__m256 y_avx = _mm256_load_ps(y);
//	//	X[i] = y;
//	//}
//
//	////float const X_out = *X;
//	//double y __attribute__((aligned(32))) = {1,2,3};
//	//__m256d sine_wave = _m256_load_pd(y);
	////_mm256_store_ps(&X[0], sine_wave);
// /	return sine_wave;
//}

__m256 generate_sine_wave(int freq, int periods, int len){
	int i,j;
	float * sin_x = new (std::align_val_t(32)) float[len];
	float * sin_y = new (std::align_val_t(32)) float[len];
	float samples = 2 * PI * periods * freq;
	
	for (i = 0; i < len; i++){
		float x = (float) (i/freq);
		sin_y[i * sizeof(float)] = sin(x);
	//	sin_x[i * sizeof(float)] = sin(x);
		//cout << a[i * sizeof(float)] << " ";
	}
	__m256 y = _mm256_load_ps(sin_y);
	
	return y;

}


int main(){
	int freq = 10;
	int periods = 1;
	int len = (int) 2 * PI * periods * freq;

	__m256 wave = generate_sine_wave(freq, periods, len);

	return 0;
}