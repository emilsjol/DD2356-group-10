#include <stdlib.h>
#include <omp.h>
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
    #pragma omp critical {
        sum_val += x[i];
    }
  }

  return sum_val;
}

double parallel_sum_3(double *x, size_t size)
{
  double sum_val = 0.0;

  double local_sum[] = malloc(MAX_THREADS * sizeof(double));
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

  double local_sum[] = malloc(SPACING*MAX_THREADS * sizeof(double));
  #pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    local_sum[omp_get_thread_num()] += x[SPACING*i];
  }

  for (size_t i = 0; i < MAX_THREADS; i++) {
    sum_val += local_sum[SPACING*i];
  }

  return sum_val;
}