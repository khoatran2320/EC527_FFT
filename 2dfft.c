#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <immintrin.h>

// gcc 2dfft.c -o 2dfft -w -mavx -lm && ./2dfft


// https://ocw.mit.edu/courses/mathematics/18-335j-introduction-to-numerical-methods-spring-2019/week-14/MIT18_335JS19_lec39.pdf
// https://www.codeproject.com/Articles/874396/Crunching-Numbers-with-AVX-and-AVX
#define		big_n 		10
#define 	RADIX		2 		
#define M_PI 3.14159265358979323846


void generate_twiddle_factors(double complex twiddle_factors[], int N, int inverse){
	for (int i = 0; i < (N/2) - 1; i++) {
		twiddle_factors[i] = cexp((inverse ? 2 : -2) * M_PI * i / N * I);
	}
}

void fft_1D(int x[], int X[], int N){
	int i,j,k;
	for (int r = 0; r < N; r++) {
		int r_even = r;
		int r_odd = r + 1;
		
	}

}


int main(){
	int x[] = {1,2,3,4};
	int r = (big_n / RADIX) - 1;
	int X[r];
	fft_1D(x, X, big_n);

	double complex *twiddle_factors = malloc((big_n) * sizeof(double complex));
	generate_twiddle_factors(twiddle_factors, big_n, 0);


	return 0;
}

