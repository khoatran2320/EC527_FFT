
#include <math.h>
#include <iostream>
#include <complex>
#include <immintrin.h>
#include <bitset>
#include <time.h>

#define PI 3.14159265358979323846
#define PERIODS 1
#define FREQ	100
using namespace std;



void transpose_matrix(float *in, unsigned long N){
	int i;
	int j;
	float temp;

	for (i = 0; i < N; i++)
	{
		for (j=i; j < N; j++)
		{
			temp = in[i*N+j];
			in[i*N+j] = in[j*N+i];
			in[j*N+i] = temp;
		}
	}
}


void init_array(float *in, int m){
	for (int i = 0; i < m; i ++){
		for (int j = 0; j < m; j++) {
			in[i * m + j] = i * m + j;
		}
	}
}
void intel_transpose(float *p){
	// https://software.intel.com/content/www/us/en/develop/articles/3d-vector-normalization-using-256-bit-intel-advanced-vector-extensions-intel-avx.html
//float *p;  // address of first vector
 __m128 *m = (__m128*) p;
 __m256 m03;
 __m256 m14; 
 __m256 m25; 
 m03  = _mm256_castps128_ps256(m[0]); // load lower halves
 m14  = _mm256_castps128_ps256(m[1]);
 m25  = _mm256_castps128_ps256(m[2]);
 m03  = _mm256_insertf128_ps(m03 ,m[3],1);  // load upper halves
 m14  = _mm256_insertf128_ps(m14 ,m[4],1);
 m25  = _mm256_insertf128_ps(m25 ,m[5],1);

 __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE( 2,1,3,2)); // upper x's and y's 
 __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE( 1,0,2,1)); // lower y's and z's
 __m256 x  = _mm256_shuffle_ps(m03, xy , _MM_SHUFFLE( 2,0,3,0)); 
 __m256 y  = _mm256_shuffle_ps(yz , xy , _MM_SHUFFLE( 3,1,2,0)); 
 __m256 z  = _mm256_shuffle_ps(yz , m25, _MM_SHUFFLE( 3,0,3,1)); 

}

double interval(struct timespec start, struct timespec end)
{
  struct timespec temp;
  temp.tv_sec = end.tv_sec - start.tv_sec;
  temp.tv_nsec = end.tv_nsec - start.tv_nsec;
  if (temp.tv_nsec < 0) {
    temp.tv_sec = temp.tv_sec - 1;
    temp.tv_nsec = temp.tv_nsec + 1000000000;
  }
  return (((double)temp.tv_sec) + ((double)temp.tv_nsec)*1.0e-9);
}

int main(){
	int m = 5;
	float arr[m * m];
//	init_array(arr, m);
	transpose_matrix(arr,m);
	intel_transpose(arr);

	return 0;
}