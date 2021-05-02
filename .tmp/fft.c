#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <immintrin.h>

#define PI 3.14159265358979323846
#define PERIODS 1
#define FREQ	100
//__m256 v=_mm256_load_ps(&f[0]);
//v=_mm256_add_ps(v,v);
//_mm256_store_ps(&f[0],v);

///https://stackoverflow.com/questions/35011579/fill-constant-floats-in-avx-intrinsics-vec
//https://www.osti.gov/servlets/purl/1514771

// gcc fft.c -o fft -w -lm -mavx2

void generate_twiddle_factors(double complex twiddle_factors[], int N, int inverse){
	for (int i = 0; i < (N/2) - 1; i++) {
		twiddle_factors[i] = cexp((inverse ? 2 : -2) * M_PI * i / N * I);
	}
}

//__m256d generate_sine_wave(float periods, float freq){
//	// X[samples] =
//
//	float X[(int) (2 * PI * periods * freq) + 1];// __attribute__((aligned(64)));	
//	__assume_aligned(X,64);
//	double samples = 2 * PI * periods * freq;
//	for (int i = 0; i < (int) samples; i++) {
//		double x = (double) (i/freq);
//		double y =  sin(x);
//
//		//__m256 y_avx = _mm256_load_ps(y);
//		X[i] = y;
//	}
//
//	//float const X_out = *X;
//	double y __attribute__((aligned(32))) = {1,2,3};
//	__m256d sine_wave = _m256_load_pd(y);
//	//_mm256_store_ps(&X[0], sine_wave);
//	return sine_wave;
//}

int main(){
	int i, j;
	int periods = 1;
	int freq = 100;
	int len = (int) (2 * PI * periods * freq) + 1;
//	float constexpr ALIGN(float a[2 * PI * periods * freq]);
	float X[(int) (2 * PI * periods * freq) + 1];// __attribute__((aligned(64)));	
	__m256 sine_wave[len];

	
	return 0;
}