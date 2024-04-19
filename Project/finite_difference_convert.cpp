#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

vector<vector<float> > matrix_add(vector<vector<float> > &matrix1, vector<vector<float> > &matrix2) {
	int rows = matrix1.size();
	int cols = matrix1[0].size();
	vector<vector<float> > temp(rows, vector<float>(cols));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			temp[i][j] = matrix1[i][j] + matrix2[i][j];
		}
	}
	return temp;
}

vector<vector<float> > matrix_scalar_multiply(vector<vector<float> > &matrix, float scalar) {
	int rows = matrix.size();
	int cols = matrix[0].size();
	vector<vector<float> > temp(rows, vector<float>(cols));
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			temp[i][j] = matrix[i][j] * scalar;
		}
	}
	return temp;
}

void roll(vector<vector<float> > &matrix, int shift_rows, int shift_cols)
{
	int rows = matrix.size();
	int cols = matrix[0].size();

	// Calculate effective row shift within range [0, rows)
	shift_rows = (shift_rows % rows + rows) % rows;

	// Calculate effective column shift within range [0, cols)
	shift_cols = (shift_cols % cols + cols) % cols;

	// Temporary matrix to hold rolled elements
	vector<vector<float> > temp(rows, vector<float>(cols));

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
int main()
{
	vector<std::vector<float> > test = {
		{0, 1, 2},
		{3, 4, 5},
		{6, 7, 8}};
	//first arg = 1 shifts vertical 1
	//first arg = -1 shifts vertical -1
	//second arg = 1 shifts horizontal 1
	//second arg = -1 shifts horizontal -1
	roll(test, 0, 0);
	vector<vector<float> > result = matrix_scalar_multiply(test, 10.0);

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
