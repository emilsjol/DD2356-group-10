#include <stdlib.h>
#include <omp.h>
#include <stdio.h>

int MAX_THREADS = 32;
void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double serial_sum(double *x, size_t size)
{
  double sum_val = 0.0;

  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
  }

  return sum_val;
}

double parallel_sum_1(double *x, size_t size)
{
  double sum_val = 0.0;

  #pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
  }

  return sum_val;
}

double parallel_sum_2(double *x, size_t size)
{
  double sum_val = 0.0;

  #pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    #pragma omp critical
    {    
	    sum_val += x[i];
    }
  }

  return sum_val;
}

double parallel_sum_3(double *x, size_t size)
{
  double sum_val = 0.0;

  double* local_sum = malloc(MAX_THREADS * sizeof(double));
  #pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    local_sum[omp_get_thread_num()] += x[i];
  }

  for (size_t i = 0; i < MAX_THREADS; i++) {
    sum_val += local_sum[i];
  }

  return sum_val;
}

double parallel_sum_4(double *x, size_t size)
{
  double sum_val = 0.0;

  int SPACING = 8;

  double* local_sum = malloc(SPACING*MAX_THREADS * sizeof(double));
  #pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    local_sum[omp_get_thread_num()*SPACING] += x[i];
  }

  for (size_t i = 0; i < MAX_THREADS; i++) {
    sum_val += local_sum[SPACING*i];
  }

  return sum_val;
}

int main(int argc, char* argv[]) {
  size_t SIZE = 100000000;
  double* input = malloc(sizeof(double) * SIZE);
  generate_random(input, SIZE);

  double sum = 0.0;
  if(*argv[1] == '0') {
    sum = serial_sum(input, SIZE);
  } else if(*argv[1] == '1') {
    sum = parallel_sum_1(input, SIZE);
  } else if(*argv[1] == '2') {
    sum = parallel_sum_2(input, SIZE);
  } else if(*argv[1] == '3') {
    sum = parallel_sum_3(input, SIZE);
  } else if(*argv[1] == '4') {
    sum = parallel_sum_4(input, SIZE);
  }

  printf("Sum is %f\n", sum);
}
