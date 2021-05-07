/*****************************************************************************/
// nvcc -arch sm_35 mmm_shared.cu -o mmm_shared

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <assert.h>
#include <cuComplex.h>

#define GIG 1000000000
#define CPG 3.0           // Cycles per GHz -- Adjust to your computer

#define TILE_WIDTH 32
#define SIZE 8192
#define BLOCK_SIZE (SIZE / 8)
#define NUM_BLOCKS2 1024
#define NUM_BLOCKS 8
#define SAMPLING_RATE 100
#define FREQUENCY 2

#define OPTIONS 1
#define IDENT 0


static void generate_sin_points_c(cuDoubleComplex *v, unsigned long size, unsigned long sampling_rate, unsigned long frequency){
    double time = 0.0;
    double inc = (double)1/sampling_rate;
    double W = (double)2 * M_PI * frequency;
    int i;

    for(i = 0; i < size; ++i){
        v[i] = make_cuDoubleComplex(100*sin(time*W), 0);
		//printf("v[%d] is %.2lf j%.2lf\n", i, creal(v[i]), cimag(v[i])); //Used for debugging, checked why vector results in function didn't match result outside function
        time += inc;
    }
}

// Assertion to check for errors
#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), (char *)__FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
  {
    fprintf(stderr, "CUDA_SAFE_CALL: %s %s %d\n",
                                       cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

//My CUDA function for global FFT that works on entire matrix
__global__ void kernel_FFT (int rowlen, cuDoubleComplex * exptable, cuDoubleComplex * fft_matrix, int num_per_thread_y) {
  
    int levels = 0;
    cuDoubleComplex temp;
    int i, j, k, l, m, size;
    int val;
    int start_of_thread = num_per_thread_y*threadIdx.y;

	for (i = rowlen; i > 1U; i >>= 1)
		levels++;
	
	// Bit-reversed addressing permutation
	int i_index;
	int j_index;
	int l_index;

    /* Determine the row number thread is acting on and precompute stride accordingly */
    int row = blockIdx.x * (rowlen / gridDim.x) + threadIdx.x;
    int stride = row * rowlen;


    /* Swap the vector elements */
    for (i = start_of_thread; i < start_of_thread + num_per_thread_y; i++) {
        val = i;
        j = 0;
        for (k = 0; k < levels; k++, val >>= 1)
            j = (j << 1) | (val & 1U);

        if (j > i) {
            i_index = stride + i;
            j_index = stride + j;
            cuDoubleComplex temp = fft_matrix[i_index];
            fft_matrix[i_index] = fft_matrix[j_index];
            fft_matrix[j_index] = temp;
        }
    }
    __syncthreads();
    
    // Cooley-Tukey decimation-in-time radix-2 FFT
    // loop through each stage

    for (size = 2; size <= rowlen; size *= 2) {

        __syncthreads();

        int halfsize = size / 2;																			
        int tablestep = rowlen / size;		
        bool skipped = 0;

        //for each stage, compute butterly for 2 outputs in groups of 2, 4, 8, ...

        if ( (halfsize >= num_per_thread_y) && ( (start_of_thread/halfsize) % 2) ) //If halfsize is larger than a thread, if the start of current thread is the second half of a butterfly, just ignore it for this iteration
        {
          continue;
        }


        for (i = start_of_thread; i < start_of_thread + num_per_thread_y; i+=size) {

            k = (i % halfsize) * tablestep; //If start of thread is in middle of a halfsize, k should be set accordingly.

            for (j = i; (j < i + halfsize && j < start_of_thread + num_per_thread_y) ; j++, k += tablestep) {   // Keep going until you either exceed the halfsize or you cross into next thread (for when halfsize is larger than thread size)
                int l = j + halfsize;

                j_index = (stride) + j;
                l_index = (stride) + l;									
                
                //printf("Row %d, thread %d: j = %d l = %d. j_index = %d, l_index = %d\n", row, threadIdx.y, j, l, j_index, l_index);
                
                temp = cuCmul( fft_matrix[l_index] , exptable[k]);
                fft_matrix[l_index] = cuCsub(fft_matrix[j_index], temp);
                fft_matrix[j_index] = cuCadd(fft_matrix[j_index], temp);
            }

            if (size == rowlen)  // Prevent overflow in 'size *= 2'
                break;
        }
    }
    __syncthreads();
}


//My CUDA function for matrix transpose (in place)
__global__ void kernel_InPlaceTranspose (int rowlen, cuDoubleComplex * transpose_matrix) {

    int col = blockIdx.x * (rowlen / gridDim.x) + threadIdx.x;
    int i;
    cuDoubleComplex temp;

    int row_index;
    int col_index;

    for (i = col; i < rowlen; i++)
    {
        row_index = (col * rowlen) + i;
        col_index = (i * rowlen) + col;
        temp = transpose_matrix[row_index];
        transpose_matrix[row_index] = transpose_matrix[col_index];
        transpose_matrix[col_index] = temp;
    }
}




//My CUDA function for matrix transpose (in place)
//Taken from https://developer.nvidia.com/blog/efficient-matrix-transpose-cuda-cc/
__global__ void kernel_Transpose (int rowlen, cuDoubleComplex * in_matrix, cuDoubleComplex * out_matrix) {

    __shared__ cuDoubleComplex tile[TILE_WIDTH][TILE_WIDTH];

    int col = blockIdx.x * TILE_WIDTH + threadIdx.x;
    int row = blockIdx.y * TILE_WIDTH + threadIdx.y;
    
    tile[threadIdx.y][threadIdx.x] = in_matrix[row*rowlen + col];
    __syncthreads();

    out_matrix[col*rowlen + row] = tile[threadIdx.y][threadIdx.x];
}


typedef cuDoubleComplex data_t;

/* Create abstract data type for matrix */
typedef struct {
  long int len;
  data_t *data;
} matrix_rec, *matrix_ptr;


/*****************************************************************************/
int main(int argc, char *argv[])
{
  int OPTION;
  double interval(struct timespec start, struct timespec end);
  struct timespec time1, time2;

  int clock_gettime(clockid_t clk_id, struct timespec *tp);
  matrix_ptr new_matrix(long int len);
  int set_matrix_row_length(matrix_ptr m, long int row_len);
  long int get_matrix_length(matrix_ptr m);
  int init_matrix(matrix_ptr m, long int len);
  int zero_matrix(matrix_ptr m, long int len);
  void mmm_kij(matrix_ptr a, matrix_ptr b, matrix_ptr c);
  double fRand(double fMin, double fMax);
  int copy_matrix(data_t *original, data_t *copy, long int MAXSIZE);
  data_t *get_matrix_start(matrix_ptr m);
  
  long int i, j, k;
  long int levels = 0;
  long int time_sec, time_ns;
  float t;

  // GPU Timing variables
  cudaEvent_t start, stop, start2, stop2;
  float elapsed_gpu, elapsed_just_FFT;

  // Arrays on GPU global memory
  data_t *FFT_gpu;
  data_t *FFT_gpu_transpose;
  data_t *FFT_host;




  /* declare and initialize the matrix structure */
  matrix_ptr fft_matrix = new_matrix(SIZE);
  fft_matrix->len = SIZE;
  generate_sin_points_c(fft_matrix->data, SIZE*SIZE, SAMPLING_RATE, FREQUENCY);

/*
    printf("\n\n\nResult of GPU code\n");  
    for(i = 0; i < SIZE; ++i){
        for (j = 0; j < SIZE; ++j){
        printf("%.2lf j%.2lf   ", cuCreal(fft_matrix->data[i*SIZE+j]), cuCimag(fft_matrix->data[i*SIZE+j]) );
        }
        printf("\n");
    }
*/


  cuDoubleComplex *exptable = (cuDoubleComplex *) malloc((SIZE / 2) * sizeof(cuDoubleComplex));
  cuDoubleComplex *exptable_gpu;

  bool inverse = 0;
  cuDoubleComplex num;
  
  //printf("\n\n\n\n");
  for (i = 0; i < SIZE / 2; i++)
  {
    num = make_cuDoubleComplex( 0, (inverse ? 2: -2) * M_PI * i / SIZE );
    t = expf(num.x);
    sincos (num.y, &exptable[i].y, &exptable[i].x);
    exptable[i].x *= t;
    exptable[i].y *= t;
    //printf("%.2lf j%.2lf     ", cuCreal(exptable[i]), cuCimag(exptable[i]) );
  }
  //printf("\n\n\n\n\n");

  for (i = SIZE; i > 1U; i >>= 1)
		levels++;

  // Allocate CPU memory
  size_t allocSize = (SIZE*SIZE) * sizeof(cuDoubleComplex);
  size_t allocSize2 = (SIZE/2) * sizeof(cuDoubleComplex);

  FFT_host          = (cuDoubleComplex *) malloc(allocSize);
  FFT_gpu           = (cuDoubleComplex *) malloc(allocSize);
  FFT_gpu_transpose = (cuDoubleComplex *) malloc(allocSize);

  // Select GPU
  CUDA_SAFE_CALL(cudaSetDevice(0));


  // Create the cuda events
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventCreate(&start2);
  cudaEventCreate(&stop2);

  // Record event on the default stream
  cudaEventRecord(start, 0);

  // Allocate GPU memory

  CUDA_SAFE_CALL(cudaMalloc((void **)&FFT_gpu, allocSize));
  CUDA_SAFE_CALL(cudaMalloc( (void **)&exptable_gpu, allocSize2 ));
  //CUDA_SAFE_CALL(cudaMalloc( (void **)&FFT_gpu_transpose, allocSize ));
  
  CUDA_SAFE_CALL(cudaMemcpy(FFT_gpu,fft_matrix->data,allocSize,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(exptable_gpu,exptable,allocSize2,cudaMemcpyHostToDevice));
  
  // Launch the kernels to make 2D FFT happen!

  int num_threads_y = 1;

  while ( (BLOCK_SIZE * (SIZE / num_threads_y)) > 1024 )
    num_threads_y <<= 1;
   
  printf("num per thread y is %d\n", num_threads_y);

  dim3 dimGrid(NUM_BLOCKS,1,1);
  dim3 dimBlock(BLOCK_SIZE,SIZE/num_threads_y,1);

  dim3 dimGrid2(NUM_BLOCKS2,1,1);
  dim3 dimBlock2(SIZE/NUM_BLOCKS2,1,1);

  cudaEventRecord(start2, 0);

  kernel_FFT <<<dimGrid, dimBlock>>>(SIZE, exptable_gpu, FFT_gpu, num_threads_y);
  CUDA_SAFE_CALL(cudaPeekAtLastError());
  
  //kernel_Transpose <<<dimGrid2, dimBlock2>>>(SIZE, FFT_gpu, FFT_gpu_transpose);
  kernel_InPlaceTranspose <<<dimGrid2, dimBlock2>>> (SIZE, FFT_gpu);
  CUDA_SAFE_CALL(cudaPeekAtLastError());
  

  kernel_FFT <<<dimGrid, dimBlock>>>(SIZE, exptable_gpu, FFT_gpu, num_threads_y);
  CUDA_SAFE_CALL(cudaPeekAtLastError());
  
  //kernel_Transpose <<<dimGrid2, dimBlock2>>>(SIZE, FFT_gpu_transpose, FFT_gpu);
  kernel_InPlaceTranspose <<<dimGrid2, dimBlock2>>> (SIZE, FFT_gpu);

  CUDA_SAFE_CALL(cudaPeekAtLastError());


  cudaEventRecord(stop2, 0);
  cudaEventSynchronize(stop2);
  cudaEventElapsedTime(&elapsed_just_FFT, start2, stop2);

  // Check for errors during launch
  CUDA_SAFE_CALL(cudaPeekAtLastError());
  CUDA_SAFE_CALL(cudaMemcpy(FFT_host,FFT_gpu,allocSize,cudaMemcpyDeviceToHost));

  //CUDA_SAFE_CALL(cudaFree(FFT_gpu_transpose));
  CUDA_SAFE_CALL(cudaFree(FFT_gpu));
  CUDA_SAFE_CALL(cudaFree(exptable_gpu));
  
  // Stop and destroy the timer
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_gpu, start, stop);
  printf("\nGPU total time: %f (msec), GPU MMM time: %f (msec)\n", elapsed_gpu, elapsed_just_FFT);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  cudaEventDestroy(start2);
  cudaEventDestroy(stop2);

    printf("[");
/*    for(i = 0; i < SIZE; ++i){
        for (j = 0; j < SIZE; ++j){
        printf("%.2lf + %.2lfi,   ", cuCreal(FFT_host[i*SIZE+j]), cuCimag(FFT_host[i*SIZE+j]) );
        }
        printf(";\n");
    }
    printf("];");
*/

}/* end main */













/**********************************************/

/* Returns a random number between fMin and fMax */
double fRand(double fMin, double fMax)
{
  double f = (double)random() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}


/* Create matrix of specified length */
matrix_ptr new_matrix(long int len)
{
  long int i;

  /* Allocate and declare header structure */
  matrix_ptr result = (matrix_ptr) malloc(sizeof(matrix_rec));
  if (!result) return NULL;  /* Couldn't allocate storage */
  result->len = len;

  /* Allocate and declare array */
  if (len > 0) {
    data_t *data = (data_t *) calloc(len*len, sizeof(cuDoubleComplex));
    if (!data) {
	  free((void *) result);
	  printf("\n COULDN'T ALLOCATE %ld BYTES STORAGE \n", result->len);
	  return NULL;  /* Couldn't allocate storage */
	}
	result->data = data;
  }
  else result->data = NULL;

  return result;
}

/* Set length of matrix */
int set_matrix_row_length(matrix_ptr m, long int row_len)
{
  m->len = row_len;
  return 1;
}

/* Return length of matrix */
long int get_matrix_length(matrix_ptr m)
{
  return m->len;
}

/* initialize matrix */
int init_matrix(matrix_ptr m, long int len)
{
  long int i;

  if (len > 0) {
    m->len = len;
    for (i = 0; i < len*len; i++)
      m->data[i] = make_cuDoubleComplex( fRand((double)(5.0),(double)(15.0)) , 0);
    return 1;
  }
  else return 0;
}

/* initialize matrix */
int zero_matrix(matrix_ptr m, long int len)
{
  long int i,j;

  if (len > 0) {
    m->len = len;
    for (i = 0; i < len*len; i++)
      m->data[i] = make_cuDoubleComplex( IDENT, 0 );
    return 1;
  }
  else return 0;
}

int copy_matrix(data_t *original, data_t *copy, long int MAXSIZE)
{
  int i;

  printf("Made it here!\n");
  for (i = 0; i < MAXSIZE * MAXSIZE; i++)
  {
    printf("%d, ", i);
    copy[i] = original[i];
  }
  printf("\n");
  return 1;
}

data_t *get_matrix_start(matrix_ptr m)
{
  return m->data;
}

/*************************************************/

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