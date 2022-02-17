#include"NKeval.hpp"

vector<vector<double> > NK_land()
{
	int num = pow(2, N);
	vector<vector<double>> binary;
	vector<double> binaryRow;
	for (int i = 0; i < num; i++) {
		//cout << i << " : ";
		for (int j = 0; j < N; j++) {
			binaryRow.push_back((double)rand() / RAND_MAX);
			//cout << binaryRow[j] << " : ";
		}
		//cout << endl;
		binary.push_back(binaryRow);
		binaryRow.clear();
	}
	return binary;
}

vector<vector<size_t>> IntMatrix()
{
	vector<vector<size_t>> Matrix;
	vector<size_t> MatrixRow;
	vector<size_t> Index1;
	mt19937 get_rand_mt(Landscape_seed);
	int temp;
	vector<size_t> choseOne;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			MatrixRow.push_back(0);
		}
		Matrix.push_back(MatrixRow);
		MatrixRow.clear();
	}
	for (int aa1 = 0; aa1 < N; aa1++) {
		for (int i = 0; i < N; i++) {
			Index1.push_back(i);
		}
		Index1.erase(Index1.begin() + aa1);
		shuffle(Index1.begin(), Index1.end(), get_rand_mt);
		Index1.push_back(aa1);
		for (int j = K; j >= 0; j--) {
			choseOne.push_back(Index1[N - j - 1]);
		}
		for (int k = 0; k <= K; k++) {
			Matrix[aa1][choseOne[k]] = 1;
		}
		Index1.clear();
		choseOne.clear();
	}
	return Matrix;
}

vector<vector<size_t>> AllBitArrays()
{
	vector<vector<size_t>> Current_position;
	vector<size_t> Current_position_Row;
	int prenum;
	int num = pow(2, N);

	//ç≈èâÇÃ0ÇÕÇ†ÇÁÇ©Ç∂ÇﬂîzóÒÇ…ì¸ÇÍÇƒÇ®Ç≠ÅB
	for (int i = 0; i < N; i++)
	{
		Current_position_Row.push_back(0);
	}
	Current_position.push_back(Current_position_Row);
	Current_position_Row.clear();
	for (int i = 1; i < num; i++)
	{
		prenum = i;
		for (int j = 0; j < N; j++)
		{
			Current_position_Row.push_back(prenum % 2);
			prenum = prenum / 2;
			//ëºÇÃbitóÒÇÕëSÇƒ0Ç…Ç∑ÇÈ
			if (prenum == 0)
			{
				for (int k = j + 1; k < N; k++)
				{
					Current_position_Row.push_back(0);
				}
				break;
			}
		}
		Current_position.push_back(Current_position_Row);
		Current_position_Row.clear();
	}

	return Current_position;
}

vector<vector<double_t>> calc_fit(vector<vector<double_t>> NK_land, vector<vector<size_t>> Inter_m, vector<vector<size_t>> Current_position)
{
	vector<vector<double_t>> Fit_vector;
	vector<double_t> Fit_Row;
	vector<size_t> sumRow;
	int number = 0;
	int num = pow(2, N);

	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				sumRow.push_back(Inter_m[j][k] * Current_position[i][k]);
				if (sumRow[k] == 1)
				{
					number += pow(2, k);
				}
			}
			Fit_Row.push_back(NK_land[number][j]);
			number = 0;
			sumRow.clear();
		}
		Fit_vector.push_back(Fit_Row);
		Fit_Row.clear();
	}
	return Fit_vector;
}

vector<double_t> comb_and_values(vector<vector<size_t>> Current_position, vector<vector<double_t>> Fit_vector)
{
	vector<double_t> evalList;
	double total = 0.0;
	int num = pow(2, N);

	for (int i = 0; i < num; i++)
	{

		for (int k = 0; k < N; k++)
		{
			total += Fit_vector[i][k];
		}
		evalList.push_back(total / (double)N);
		total = 0.0;
	}
	return evalList;
}
