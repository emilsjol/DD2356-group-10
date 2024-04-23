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

vector<vector<float>  > laplacian(vector<vector<float>  > &ulx, vector<vector<float>  > &uly,
								vector<vector<float>  > &urx, vector<vector<float>  > &ury, vector<vector<float>  > &u)
{
	vector<vector<float>  > u4 = matrix_scalar_multiply(u, -4.0);

	vector<vector<float>  > uX = matrix_add(ulx, urx);

	vector<vector<float>  > uY = matrix_add(uly, ury);

	vector<vector<float>  > uRes = matrix_add(uX, uY);

	vector<vector<float>  > result = matrix_add(u4, uRes);

	return result;
}

void roll(vector<vector<float>  > &matrix, int shift_rows, int shift_cols)
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

	matrix = temp;
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
	int n = 256; //resolution
	int boxsize = 1;
	int c = 1;
	int t = 0;
	int t_end = 2;
	int plot_real_time = true;

	float dx = boxsize / n;
	float dt = (sqrt(2.0) / 2.0 ) * dx / c;
	int aX = 0;
	int aY = 1;
	int r = -1;
	int l = 1;
	float fac = (dt*dt) * (c*c) / (dx*dx);
	
	vector<float> linspace = create_lin_space(0.5*dx, boxsize-(0.5*dx), n);
	//Y, X = ??

	vector<vector<float> > U = create_zero_matrix(n);
	//go crazy with masks
	vector<vector<bool> > mask = create_bool_matrix(n);
	
	for (int i = 0; i < n; i++) {
		mask[0][i] = true;
	}
	for (int i = 0; i < n; i++) {
		mask[n-1][i] = true;
	}
	for (int i = 0; i < n; i++) {
		mask[i][0] = true;
	}
	for (int i = 0; i < n; i++) {
		mask[n-1][0] = true;
	}
	for (int i = 0; i < n; i++) {
		for (int j = int(n/4); j < int(n*9/32); j++) {
			mask[j][i] = true;
		}
	}
	for (int i = int(n*5/16); i < int(n*3/8); i++) {
		for (int j = 1; j < n; j++) {
			mask[j][i] = false;
		}
	}
	for (int i = int(n*5/8); i < int(n*11/16); i++) {
		for (int j = 1; j < n; j++) {
			mask[j][i] = false;
		}
	}

	vector<vector<float> > Uprev = matrix_scalar_multiply(U, 1.0);



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

	result = laplacian(test, test, test, test, test);
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
