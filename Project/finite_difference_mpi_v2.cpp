#include <cmath>
#include <vector>
#include <iostream>
#include <math.h>
//#include <omp.h>
#include <mpi.h>
using namespace std;

int num_processes = 0;

vector<double> matrix_add(vector<double> &matrix1, vector<double> &matrix2)
{
	int length = matrix1.size();
    vector<double> matrix1_recv(length / num_processes);
    MPI_Scatter(&matrix1[0], length / num_processes, MPI_DOUBLE, &matrix1_recv[0], length / num_processes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    vector<double> matrix2_recv(length / num_processes);
    MPI_Scatter(&matrix2[0], length / num_processes, MPI_DOUBLE, &matrix2_recv[0], length / num_processes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	vector<double> temp(length);
	for (int i = 0; i < length; i++)
	{
        temp[i] = matrix1_recv[i] + matrix2_recv[i];
	}
    vector<double> temp_recv(length); 
    MPI_Gather(&temp[0], length / num_processes, MPI_DOUBLE, &temp_recv[0], length / num_processes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return temp_recv;
}

vector<double> matrix_scalar_multiply(vector<double> &matrix, double scalar)
{
    vector<double> matrix_recv(matrix.size() / num_processes);
    MPI_Scatter(&matrix[0], matrix.size() / num_processes, MPI_DOUBLE, &matrix_recv[0], matrix.size() / num_processes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    int length = matrix_recv.size();
	vector<double> temp(length);
	for (int i = 0; i < length; i++)
	{
        temp[i] = matrix_recv[i] * scalar;
	}
    vector<double> temp_recv(matrix.size()); 
    MPI_Gather(&temp[0], matrix.size() / num_processes, MPI_DOUBLE, &temp_recv[0], matrix.size() / num_processes, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	return temp_recv;
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

vector<double> roll(vector<double> &matrix, int shift_rows, int shift_cols, int rows, int cols, int rank)
{
    int length = matrix.size();
    int per_rank = length / num_processes;

    // Different logic due to ghost rows required/not

    if(shift_rows == 1) {
        // Shift rows ==> Required ghost row
        vector<double> matrix_recv(per_rank + shift_rows*(cols));

        // Receive into a different index, no shifting operations needed!
        MPI_Scatter(&matrix[0], per_rank, MPI_DOUBLE, &matrix_recv[cols], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(rank == 0) 
        {
            // Send ghost row from rank 0. Serial, but slightly easier to implement this way
            for(int i = 1; i < num_processes; i++) 
            {
                MPI_Ssend(&matrix[(i * per_rank) - cols], cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            // Rank 0 can just read instead of sending to itself
            for(int i = 0; i < cols; i++) {
                matrix_recv[i] = matrix[(num_processes * per_rank) - cols + i];
            }
        } 
        else
        {
            MPI_Status status;
            // Receive into empty part of receiving matrix
            MPI_Recv(&matrix_recv[0], cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
        vector<double> temp_recv(rows*cols);
        MPI_Gather(&matrix_recv[0], per_rank, MPI_DOUBLE, &temp_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return temp_recv;
    }

    if(shift_rows == -1) {
        // Shift rows ==> Required ghost row
        vector<double> matrix_recv(per_rank + cols);

        // Receive into index 0, but allocate more memory after for the ghost row
        MPI_Scatter(&matrix[0], per_rank, MPI_DOUBLE, &matrix_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if(rank == 0) 
        {
            // Send ghost rows serially from rank 0
            for(int i = 1; i < num_processes; i++) 
            {
                MPI_Ssend(&matrix[(((i+1) % num_processes) * per_rank)], cols, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            // Rank 0 fixes for itself
            for(int i = 0; i < cols; i++) {
                matrix_recv[per_rank + i] = matrix[per_rank + i];
            }
        } 
        else
        {
            MPI_Status status;
            // Receive ghost row
            MPI_Recv(&matrix_recv[per_rank], cols, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }
        vector<double> temp_recv(rows*cols);
        // Sends from a shifted pointer, thereby no shifting required
        MPI_Gather(&matrix_recv[cols], per_rank, MPI_DOUBLE, &temp_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return temp_recv;
    }

    if (shift_cols == -1) {
        // Column shift. No ghost cells required
        vector<double> matrix_recv(per_rank+1);
        MPI_Scatter(&matrix[0], per_rank, MPI_DOUBLE, &matrix_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Move the edge column
        for(int i = per_rank/cols; i >= 0; i--) 
        {
            matrix_recv[(i*cols)] = matrix_recv[(i-1)*cols];
        }

        vector<double> temp_recv(rows*cols);
        // Shift everything by 1 by pointer arithmetic, fixes the remaining columns
        MPI_Gather(&matrix_recv[1], per_rank, MPI_DOUBLE, &temp_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return temp_recv;
    }

    if (shift_cols == 1) {
    
        vector<double> matrix_recv(per_rank+1);
        // Receive into an offset pointer, fixing all but one edge column
        MPI_Scatter(&matrix[0], per_rank, MPI_DOUBLE, &matrix_recv[1], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Fix the edge column
        for(int i = 0; i < per_rank / cols; i++) 
        {
            matrix_recv[i*cols] = matrix_recv[(i+1)*cols];
        }

        vector<double> temp_recv(rows*cols);
        MPI_Gather(&matrix_recv[0], per_rank, MPI_DOUBLE, &temp_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        return temp_recv;
    }

    
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
	for (int i = 0; i < (int)(n); i++) {
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

    num_processes = size;

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
	vector<double> ULX = vector<double>(N*N);
	vector<double> URX = vector<double>(N*N);
	vector<double> ULY = vector<double>(N*N);
	vector<double> URY = vector<double>(N*N);
	int counter = 0;
	while (t < tEnd) {
		counter++;
        ULX = roll(U, L, 0, N, N, rank);
        URX = roll(U, R, 0, N, N, rank);
        ULY = roll(U, 0, L, N, N, rank);
        URY = roll(U, 0, R, N, N, rank);

		vector<double> laplacian = create_laplacian(ULX, ULY, URX, URY, U);
		vector<double> twoU = matrix_scalar_multiply(U, 2.0); //2*U
		vector<double> negativeUprev = matrix_scalar_multiply(Uprev, -1.0);
		vector<double> facLaplacian = matrix_scalar_multiply(laplacian, fac);
		vector<double> first_operation = matrix_add(twoU, negativeUprev); //calculate 2*U - Uprev
		vector<double> Unew = matrix_add(first_operation, facLaplacian); //calculate 2*U - Uprev + (fac*laplacian)
		Uprev = matrix_scalar_multiply(U, 1.0);
		U = matrix_scalar_multiply(Unew, 1.0);

        int per_rank = N*N/size;	
		
        vector<double> U_recv(per_rank);

        MPI_Scatter(&U[0], per_rank, MPI_DOUBLE, &U_recv[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (int i = 0; i < per_rank; i++) {
            if (mask[(rank * (per_rank)) + i]) {
                U_recv[i] = 0.0;
            }
        }
        MPI_Gather(&U_recv[0], per_rank, MPI_DOUBLE, &U[0], per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
        

        		
        vector<double> U_recv2(per_rank);
		for (int i = 0; i < per_rank/N; i++) {
			U_recv2[i] = sin(20*M_PI*t) * (sin(M_PI*xlin[(rank*per_rank/N) + i]) * sin(M_PI*xlin[(rank*per_rank/N) + i]));
		}
        MPI_Gather(&U_recv2[0], per_rank/N, MPI_DOUBLE, &U[0], per_rank/N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  

        MPI_Bcast(&U[0], N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		t += dt;
		if(rank == 0) {
			cout << t << "\n";
			//print_matrix(U, N);
		}

	}

    if(rank == 0) {
        print_matrix(U, N);
    }
    
    MPI_Finalize();

	return 0;
}
