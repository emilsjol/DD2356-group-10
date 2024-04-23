#include <cmath>
#include <vector>
#include <iostream>
#include <math.h>
using namespace std;

vector<vector<float>  > matrix_add(vector<vector<float>  > &matrix1, vector<vector<float>  > &matrix2)
{
	int rows = matrix1.size();
	int cols = matrix1[0].size();
	vector<vector<float>  > temp(rows, vector<float>(cols));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			temp[i][j] = matrix1[i][j] + matrix2[i][j];
		}
	}
	return temp;
}

vector<vector<float>  > matrix_scalar_multiply(vector<vector<float>  > &matrix, float scalar)
{
	int rows = matrix.size();
	int cols = matrix[0].size();
	vector<vector<float>  > temp(rows, vector<float>(cols));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			temp[i][j] = matrix[i][j] * scalar;
		}
	}
	return temp;
}

vector<vector<float>  > create_laplacian(vector<vector<float>  > &ulx, vector<vector<float>  > &uly,
								vector<vector<float>  > &urx, vector<vector<float>  > &ury, vector<vector<float>  > &u)
{
	vector<vector<float>  > u4 = matrix_scalar_multiply(u, -4.0);

	vector<vector<float>  > uX = matrix_add(ulx, urx);

	vector<vector<float>  > uY = matrix_add(uly, ury);

	vector<vector<float>  > uRes = matrix_add(uX, uY);

	vector<vector<float>  > result = matrix_add(u4, uRes);

	return result;
}

vector<vector<float> > roll(vector<vector<float>  > &matrix, int shift_rows, int shift_cols)
{
	int rows = matrix.size();
	int cols = matrix[0].size();

	// Calculate effective row shift within range [0, rows)
	shift_rows = (shift_rows % rows + rows) % rows;

	// Calculate effective column shift within range [0, cols)
	shift_cols = (shift_cols % cols + cols) % cols;

	// Temporary matrix to hold rolled elements
	vector<vector<float>  > temp(rows, vector<float>(cols));

	// Roll rows
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			temp[(i + shift_rows) % rows][(j + shift_cols) % cols] = matrix[i][j];
		}
	}

	return temp;
}

vector<vector<float> > create_zero_matrix(int n) 
{
	vector<vector<float> > zero_matrix(n, vector<float>(n));
	return zero_matrix;
}

vector<vector<bool> > create_bool_matrix(int n)
{
	vector<vector<bool> > bool_matrix(n, vector<bool>(n));
	return bool_matrix;
}

vector<float> create_lin_space(float start, float end, float n)
{
	float distance = end - start;
	float increment = distance / n;
	vector<float> lin_space(n);
	for (int i = 0; i < n; i++) {
		lin_space[i] = start + (increment * i);
	}
	return lin_space;
}

int main()
{

	//start main
	int N = 256; //resolution
	int boxsize = 1;
	int c = 1;
	float t = 0;
	float tEnd = 2;
	int plotRealTime = true;

	float dx = float(boxsize) / float(N);
	float dt = (sqrt(2.0) / 2.0 ) * dx / float(c);
	cout << dt << " ";
	int aX = 0;
	int aY = 1;
	int R = -1;
	int L = 1;
	float fac = (dt*dt) * (c*c) / (dx*dx);
	
	vector<float> xlin = create_lin_space(0.5*dx, boxsize-(0.5*dx), N);
	//Y, X = ??

	vector<vector<float> > U = create_zero_matrix(N);
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
		mask[N-1][0] = true;
	}
	for (int i = 0; i < N; i++) {
		for (int j = int(N/4); j < int(N*9/32); j++) {
			mask[j][i] = true;
		}
	}
	for (int i = int(N*5/16); i < int(N*3/8); i++) {
		for (int j = 1; j < N; j++) {
			mask[j][i] = false;
		}
	}
	for (int i = int(N*5/8); i < int(N*11/16); i++) {
		for (int j = 1; j < N; j++) {
			mask[j][i] = false;
		}
	}

	vector<vector<float> > Uprev = matrix_scalar_multiply(U, 1.0);

	//main loop time
	while (t < tEnd) {
		vector<vector<float> > ULX = roll(U, L, aX);
		vector<vector<float> > URX = roll(U, R, aX);
		vector<vector<float> > ULY = roll(U, L, aY);
		vector<vector<float> > URY = roll(U, R, aY);
		vector<vector<float> > laplacian = create_laplacian(ULX, ULY, URX, URY, U);
		//update U
		vector<vector<float> > twoU = matrix_scalar_multiply(U, 2.0); //2*U
		vector<vector<float> > facLaplacian = matrix_scalar_multiply(laplacian, fac);
		vector<vector<float> > negativeUprev = matrix_scalar_multiply(Uprev, -1.0);

		vector<vector<float> > first_operation = matrix_add(twoU, negativeUprev); //calculate 2*U - Uprev
		vector<vector<float> > Unew = matrix_add(first_operation, facLaplacian); //calculate 2*U - Uprev + (fac*laplacian)
		Uprev = matrix_scalar_multiply(U, 1.0);
		U = matrix_scalar_multiply(Unew, 1.0);	
		
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (mask[j][i] == true) {
					U[j][i] = 0.0;
				}
			}
		}
		
		for (int i = 0; i < N; i++) {
			U[0][i] = sin(20*M_PI*t) * (sin(M_PI*xlin[i]) * sin(M_PI*xlin[i]));
		}
		cout << t << " ";
		t += dt;
	}


	vector<std::vector<float>  > test = {
		{0, 1, 2},
		{3, 4, 5},
		{6, 7, 8}};
	// first arg = 1 shifts vertical 1
	// first arg = -1 shifts vertical -1
	// second arg = 1 shifts horizontal 1
	// second arg = -1 shifts horizontal -1
	roll(test, 0, 0);
	vector<vector<float>  > result = matrix_scalar_multiply(test, 10.0);

	for (const auto &row : result)
	{
		for (float num : row)
		{
			cout << num << " ";
		}
		cout << endl;
	}

	result = matrix_add(test, test);

	for (const auto &row : result)
	{
		for (float num : row)
		{
			cout << num << " ";
		}
		cout << endl;
	}

	result = create_laplacian(test, test, test, test, test);
	for (const auto &row : result) 
	{
		for (float num : row)
		{
			cout << num << " ";
		}
		cout << endl;
	}
	return 0;
}
