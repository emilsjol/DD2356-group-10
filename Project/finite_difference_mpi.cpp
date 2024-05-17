#include <cmath>
#include <vector>
#include <iostream>
#include <math.h>
#include <omp.h>
#include <mpi.h>
using namespace std;

vector<vector<double>  > matrix_add(vector<vector<double>  > &matrix1, vector<vector<double>  > &matrix2)
{
	int rows = matrix1.size();
	int cols = matrix1[0].size();
	vector<vector<double>  > temp(rows, vector<double>(cols));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			temp[i][j] = matrix1[i][j] + matrix2[i][j];
		}
	}
	return temp;
}

vector<vector<double>  > matrix_scalar_multiply(vector<vector<double>  > &matrix, double scalar)
{
	int rows = matrix.size();
	int cols = matrix[0].size();
	vector<vector<double>  > temp(rows, vector<double>(cols));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			temp[i][j] = matrix[i][j] * scalar;
		}
	}
	return temp;
}
		//vector<vector<double> > laplacian = create_laplacian(ULX, ULY, URX, URY, U);
vector<vector<double>  > create_laplacian(vector<vector<double>  > &ulx, vector<vector<double>  > &uly,
								vector<vector<double>  > &urx, vector<vector<double>  > &ury, vector<vector<double>  > &u)
{
	vector<vector<double>  > u4 = matrix_scalar_multiply(u, -4.0);

	vector<vector<double>  > uX = matrix_add(ulx, urx);

	vector<vector<double>  > uY = matrix_add(uly, ury);

	vector<vector<double>  > uRes = matrix_add(uX, uY);

	vector<vector<double>  > result = matrix_add(u4, uRes);

	return result;
}

vector<vector<double> > roll(vector<vector<double>  > &matrix, int shift_rows, int shift_cols)
{
	int rows = matrix.size();
	int cols = matrix[0].size();

	// Calculate effective row shift within range [0, rows)
	shift_rows = (shift_rows % rows + rows) % rows;

	// Calculate effective column shift within range [0, cols)
	shift_cols = (shift_cols % cols + cols) % cols;

	// Temporary matrix to hold rolled elements
	vector<vector<double>  > temp(rows, vector<double>(cols));

	// Roll rows
	#pragma omp parallel for
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			temp[(i + shift_rows) % rows][(j + shift_cols) % cols] = matrix[i][j];
		}
	}

	return temp;
}

vector<vector<double> > create_zero_matrix(int n) 
{
	vector<vector<double> > zero_matrix(n, vector<double>(n));
	return zero_matrix;
}

vector<vector<bool> > create_bool_matrix(int n)
{
	vector<vector<bool> > bool_matrix(n, vector<bool>(n));
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

void print_matrix(vector<vector<double> > matrix) {
	int size = matrix.size();
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << "\n";
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
	int N = 256; //resolution
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

	vector<vector<double> > U = create_zero_matrix(N);
	//go crazy with masks
	vector<vector<bool> > mask = create_bool_matrix(N);
	
	for (int i = 0; i < N; i++) {
		mask[0][i] = true;
	}
	for (int i = 0; i < N; i++) {
		mask[N-1][i] = true;
	}
	for (int i = 0; i < N; i++) {
		mask[i][0] = true;
	}
	for (int i = 0; i < N; i++) {
		mask[i][N-1] = true;
	}
	for (int i = 0; i < N-1; i++) {
		for (int j = int(double(N)/4.0); j < int(double(N)*9.0/32.0); j++) {
			mask[j][i] = true;
		}
	}
	for (int i = int(double(N)*5.0/16.0); i < int(double(N)*3.0/8.0); i++) {
		for (int j = 1; j < N-1; j++) {
			mask[j][i] = false;
		}
	}
	for (int i = int(double(N)*5.0/8.0); i < int(double(N)*11.0/16.0); i++) {
		for (int j = 1; j < N-1; j++) {
			mask[j][i] = false;
		}
	}

	vector<vector<double> > Uprev = matrix_scalar_multiply(U, 1.0);

	//main loop time
	//MPI_Recv(&f[0], 1, MPI_DOUBLE, prev_rank, 0, MPI_COMM_WORLD, &status);
	MPI_Request * requests = new MPI_Request[size - 1];
	MPI_Request request;
	MPI_Status status;
	MPI_Status * statuses = new MPI_Status[size - 1];
	double blabla = 0.2;
	vector<vector<double> > ULX = vector<vector<double> > (N, vector<double>(N));
	vector<vector<double> > URX = vector<vector<double> > (N, vector<double>(N));
	vector<vector<double> > ULY = vector<vector<double> > (N, vector<double>(N));
	vector<vector<double> > URY = vector<vector<double> > (N, vector<double>(N));
	
	while (t < tEnd) {

		if (rank == 0) {
			cout << "Entered rank 0\n";
			ULX = roll(U, L, 0);
			cout << "Rolled rank 0\n";
		}
		if (rank == 1) {
			cout << "Entered rank 1\n";
			URX = roll(U, R, 0);
			cout << "Rolled rank 1\n";
			//MPI_Isend(&URX, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests[0]);
			MPI_Ssend(&URX, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			cout << "Semt rank 1\n";
		}
		if (rank == 2) {
			cout << "Entered rank 2\n";
			ULY = roll(U, 0, L);
			cout << "Rolled rank 2\n";
			//MPI_Isend(&ULY, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests[1]);
			MPI_Ssend(&ULY, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			cout << "Semt rank 2\n";
		}
		if (rank == 3) {
			cout << "Entered rank 3\n";
			URY = roll(U, 0, R);
			cout << "Rolled rank 3\n";
			//MPI_Isend(&URY, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests[2]);
			MPI_Ssend(&URY, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			cout << "Semt rank 3\n";
		}
		if (rank == 0) {
			cout << "Entered rank 0 receive\n";
			MPI_Recv(&blabla, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&blabla, 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&blabla, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &status);
			cout << "Exited rank 0 receive\n";
		}
		/*
		if (rank == 0) {
			cout << "Entered rank 0 receive\n";
			MPI_Irecv(&URX, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &requests[0]);
			MPI_Irecv(&ULY, 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, &requests[1]);
			MPI_Irecv(&URY, 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, &requests[2]);
		}
		*/
		//wait here for above matrices
		//MPI_Waitall(6, requests, statuses);
		cout << "BANANAS: " << rank << "\n";
		//börja parallelblock här
		vector<vector<double> > laplacian = create_laplacian(ULX, ULY, URX, URY, U);
		vector<vector<double> > twoU = matrix_scalar_multiply(U, 2.0); //2*U
		vector<vector<double> > negativeUprev = matrix_scalar_multiply(Uprev, -1.0);
		//vänta in parallel här
		//börja parallelblock här
		vector<vector<double> > facLaplacian = matrix_scalar_multiply(laplacian, fac);
		vector<vector<double> > first_operation = matrix_add(twoU, negativeUprev); //calculate 2*U - Uprev
		//vänta parallellblock här
		vector<vector<double> > Unew = matrix_add(first_operation, facLaplacian); //calculate 2*U - Uprev + (fac*laplacian)
		Uprev = matrix_scalar_multiply(U, 1.0);
		U = matrix_scalar_multiply(Unew, 1.0);	
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (mask[i][j]) {
					U[i][j] = 0.0;
				}
			}
		}
		
		//parallelisera potentiellt? vi kan kolla prestanda kanske
		for (int i = 0; i < N; i++) {
			U[0][i] = sin(20*M_PI*t) * (sin(M_PI*xlin[i]) * sin(M_PI*xlin[i]));
		}

		t += dt;
		cout << t << "\n";

	}
	return 0;
}
