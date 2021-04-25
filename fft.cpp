#include <math.h>
#include <iostream>
#include <complex>
#include <immintrin.h>

#define PI 3.14159265358979323846
#define PERIODS 1
#define FREQ	100
using namespace std;
// g++ fft.cpp -o fft_cpp -w -lm -mavx -std=c++17 -m64 && ./fft_cpp
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

// https://stackoverflow.com/questions/45554873/using-values-from-m256i-to-access-an-array-efficiently-simd/45555691
// https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/intrinsics/intrinsics-for-intel-advanced-vector-extensions/intrinsics-for-miscellaneous-operations-3/mm256-set-epi8-16-32-64x.html

void generate_sine_wave(float y[], int freq, int periods, int len){
	int i,j;
	float * sin_x = new (std::align_val_t(32)) float[len];
	float * sin_y = new (std::align_val_t(32)) float[len];
	float samples = 2 * PI * periods * freq;
	
	
	for (i = 0; i < len; i++){
		float x = ((float) i)/freq;
//		cout << x << " : " << sin(x) << endl;
//		sin_y[i * sizeof(float)] = sin(x);
		//sin_y[i] = sin(x);
		const float res = sin(x);
		y[i] = res;
		//_mm256_store_ps()
		//y[i * sizeof(float)] = _mm256_load_ps(&res);
		//cout << a[i * sizeof(float)] << " ";
		//cout << sin_y[i] << endl;
	}
	//__m256 y = _mm256_load_ps(sin_y);

	//_mm256_store_ps(sin_y, y);
//	return y;

}

static size_t reverse_bits(size_t val, int width) {
	size_t result = 0;
	int i;
	for (i = 0; i < width; i++, val >>= 1)
		result =  (result << 1 ) | (val & 1U);
	return result;

}
__m256 avx_fft(__m256 inp, size_t n, bool inverse, __m256 exptable_re, __m256 exptable_im){
	size_t temp, i, j, size, k, m;

	__m256 X;
	return X;

}
int main(){
	int freq = 10;
	int periods = 1;
	int len = (int) 2 * PI * periods * freq;
	//__m256 y[len];
	float y[len];
	//__m256 wave = generate_sine_wave(y,freq, periods, len);

	generate_sine_wave(y,freq, periods, len);

	float * fWave = new (std::align_val_t(32)) float[len];

	__m256 entry;

	//float ytest[3] = {0};
	float * ytest = new (std::align_val_t(32)) float[3];	



//	_mm256_store_ps(fWave, y);
	for (int i = 0; i < 3; i++) 
			{
				//cout << ((float)i)/freq << " : " << fWave[i] << endl;

		
	}
	cout << endl;
	int indices[8] = {0,1,2,3,4,5,6,7};

	//__m256i ind = _mm256_set_epi32(0,1,2,3,4,5,6,7);
	//__m256 result = _mm256_i32gather_ps(wave, ind, 4);

	//__mm256 index = 
	size_t n = len * len;
	int levels = 0;
	size_t temp, i, j, size, k, m;

	for (temp = n; temp > 1U; temp >>=1){
		levels++;
	}

	if ((size_t)1U << levels != n) 
		false;

	size_t i_index, j_index, l_index;

	for (m = 0; m < n; m++){
		for (i = 0; i < n; i++) {
			size_t j = reverse_bits(i, levels);

			if (j > 1) {
				i_index = m * n + i;
				j_index = m * n + j;
				//__m256 temp = __mm256_load_ps()
			}
		}
	}


	return 0;
}