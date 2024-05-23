#include <cmath>
#include <vector>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <mpi.h>
using namespace std;

vector<double> matrix_add(vector<double> &matrix1, vector<double> &matrix2)
{
	int length = matrix1.size();
	vector<double> temp(length);
	for (int i = 0; i < length; i++)
	{
        temp[i] = matrix1[i] + matrix2[i];
	}
	return temp;
}

vector<double> matrix_scalar_multiply(vector<double> &matrix, double scalar)
{
    int length = matrix.size();
	vector<double> temp(length);
	for (int i = 0; i < length; i++)
	{
        temp[i] = matrix[i] * scalar;
	}
	return temp;
}

vector<double> create_laplacian(vector<double> &ulx, vector<double> &uly,
								vector<double> &urx, vector<double> &ury, vector<double> &u)
{
	vector<double> u4 = matrix_scalar_multiply(u, -4.0);

	vector<double> uX = matrix_add(ulx, urx);

	vector<double> uY = matrix_add(uly, ury);

	vector<double> uRes = matrix_add(uX, uY);

	vector<double> result = matrix_add(u4, uRes);

	return result;
}

int index_from_two(int i, int j, int cols) {
    return ((i * cols) + j);
}

vector<double> roll(vector<double> &matrix, int shift_rows, int shift_cols, int rows, int cols)
{
    int length = matrix.size();
	// Calculate effective row shift within range [0, rows)
	shift_rows = (shift_rows % rows + rows) % rows;

	// Calculate effective column shift within range [0, cols)
	shift_cols = (shift_cols % cols + cols) % cols;

	// Temporary matrix to hold rolled elements
	vector<double> temp(length);

	// Roll rows
	for (int i = 0; i < length; ++i)
	{
        //temp[(i + shift_rows) % rows][(j + shift_cols) % cols] = matrix[i][j];
        //temp[indexFromTwo((i + shift_rows) % rows, (j + shift_cols) % cols, cols)] = matrix[i];
        temp[((((i / cols) + shift_rows) % rows) * cols) + (((i % cols) + shift_cols) % cols)] = matrix[i];
	}

	return temp;
}

vector<double> create_zero_matrix(int n) 
{
	vector<double> zero_matrix(n*n);
	return zero_matrix;
}

vector<bool> create_bool_matrix(int n)
{
	vector<bool> bool_matrix(n*n);
	return bool_matrix;
}

vector<double> create_lin_space(double start, double end, double n)
{
	double distance = end - start;
	double increment = distance / (n - 1);
	vector<double> lin_space(n);
	for (int i = 0; i < (int)(n); i++) { // Casting due to OpenMP Canonical Loop Form or something
		lin_space[i] = start + (increment * i);
	}
	return lin_space;
}

void print_matrix(vector<double> matrix, int cols) {
	cout << "Printing matrix\n";
    int length = matrix.size();
    for (int i = 0; i < length; i++) {
        cout << matrix[i] << " ";
        
        if(((i+1) % cols) == 0) 
        {
            cout << "\n";
        }
    }
    cout << "------------\n";
	return;
}

void print_vector(vector<double> v) {
	cout << "Printing vector\n";
		int size = v.size();
		for (int i = 0; i < size; i++) {
			cout << v[i] << " ";
		}
		cout << "\n------------\n";
	return;
}

int main(int argc, char* argv[])
{
	int size, rank;
	int provided = 0;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//start main
	int N = 8; //resolution
	int boxsize = 1;
	int c = 1;
	double t = 0;
	double tEnd = 2;
	int plotRealTime = true;

	double dx = double(boxsize) / double(N);
	double dt = (sqrt(2.0) / 2.0 ) * dx / double(c);

	int aX = 0;
	int aY = 1;
	int R = -1;
	int L = 1;
	double fac = (dt*dt) * (c*c) / (dx*dx);
	
	vector<double> xlin = create_lin_space(0.5*dx, boxsize-(0.5*dx), N);

	vector<double> U = create_zero_matrix(N);
	//go crazy with masks
	vector<bool> mask = create_bool_matrix(N);
	
	for (int i = 0; i < N; i++) {
		mask[index_from_two(0, i, N)] = true;
	}
	for (int i = 0; i < N; i++) {
		mask[index_from_two(N-1, i, N)] = true;
	}
	for (int i = 0; i < N; i++) {
		mask[index_from_two(i, 0, N)] = true;
	}
	for (int i = 0; i < N; i++) {
		mask[index_from_two(i, N-1, N)] = true;
	}
	for (int i = 0; i < N-1; i++) {
		for (int j = int(double(N)/4.0); j < int(double(N)*9.0/32.0); j++) {
			mask[index_from_two(j, i, N)] = true;
		}
	}
	for (int i = int(double(N)*5.0/16.0); i < int(double(N)*3.0/8.0); i++) {
		for (int j = 1; j < N-1; j++) {
			mask[index_from_two(j, i, N)] = false;
		}
	}
	for (int i = int(double(N)*5.0/8.0); i < int(double(N)*11.0/16.0); i++) {
		for (int j = 1; j < N-1; j++) {
			mask[index_from_two(j, i, N)] = false;
		}
	}

	vector<double> Uprev = matrix_scalar_multiply(U, 1.0);

	//main loop time
	MPI_Request * requests = new MPI_Request[size - 1];
	MPI_Request request;
	MPI_Status status;
	MPI_Status * statuses = new MPI_Status[size - 1];
	double blabla = 0.2;
	vector<double> ULX = vector<double>(N*N);
	vector<double> URX = vector<double>(N*N);
	vector<double> ULY = vector<double>(N*N);
	vector<double> URY = vector<double>(N*N);
	int counter = 0;
	while (t < tEnd) {
		counter++;
		if (rank == 0) {
			ULX = roll(U, L, 0, N, N);
		}
		if (rank == 1) {
			URX = roll(U, R, 0, N, N);
			MPI_Ssend(&URX[0], N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		if (rank == 2) {
			ULY = roll(U, 0, L, N, N);
			MPI_Ssend(&ULY[0], N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		if (rank == 3) {
			URY = roll(U, 0, R, N, N);
			MPI_Ssend(&URY[0], N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
		if (rank == 0) {
            MPI_Recv(&URX[0], N*N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&ULY[0], N*N, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&URY[0], N*N, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &status);
		}

		//börja parallelblock här
		vector<double> laplacian = create_laplacian(ULX, ULY, URX, URY, U);
		vector<double> twoU = matrix_scalar_multiply(U, 2.0); //2*U
		vector<double> negativeUprev = matrix_scalar_multiply(Uprev, -1.0);
		//vänta in parallel här
		//börja parallelblock här
		vector<double> facLaplacian = matrix_scalar_multiply(laplacian, fac);
		vector<double> first_operation = matrix_add(twoU, negativeUprev); //calculate 2*U - Uprev
		//vänta parallellblock här
		vector<double> Unew = matrix_add(first_operation, facLaplacian); //calculate 2*U - Uprev + (fac*laplacian)
		Uprev = matrix_scalar_multiply(U, 1.0);
		U = matrix_scalar_multiply(Unew, 1.0);	
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (mask[index_from_two(i, j, N)]) {
					U[index_from_two(i, j, N)] = 0.0;
				}
			}
		}
		
		//parallelisera potentiellt? vi kan kolla prestanda kanske
		for (int i = 0; i < N; i++) {
			U[index_from_two(0, i, N)] = sin(20*M_PI*t) * (sin(M_PI*xlin[i]) * sin(M_PI*xlin[i]));
		}

		if (rank == 0) {
				for(int i = 1; i < size; i++) {
					MPI_Ssend(&U[0], N*N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				}
		} else {
			MPI_Recv(&U[0], N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}

		t += dt;
		if(rank == 0) {
			cout << t << "\n";
			//print_matrix(U, N);
		}

	}

    /*if(rank == 0) 
    {
        cout << "Printing mask\n";
        int length = mask.size();
        for (int i = 0; i < length; i++) {
            cout << (mask[i]? "1" : "0") << " ";
            
            if(((i+1) % N) == 0) 
            {
                cout << "\n";
            }
        }
        cout << "------------\n";
    }*/

	return 0;
}
