#include <stdio.h>
#include <omp.h>

int main(int argc, char* argv[]) {

	//omp_set_num_threads(4);

	#pragma omp parallel num_threads(4)
	{
		int id = omp_get_thread_num();
		int np = omp_get_num_threads();

		printf("Hello World %d of %d\n", id, np);
	}

	return 0;
}
