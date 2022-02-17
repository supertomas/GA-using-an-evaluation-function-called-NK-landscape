#include"GAalgorithm.hpp"
#include"NKeval.hpp"

vector<vector<size_t>> init()
{
	int i, j;
	vector<vector<size_t>> genom;
	vector<size_t> genomRow; // 行用

	for (i = 0; i < MAX_GENOM_LIST; i++)
	{
		for (j = 0; j < N; j++)
		{
			if ((double)rand() / RAND_MAX <= 0.5)
			{
				genomRow.push_back(0);
			}
			else
			{
				genomRow.push_back(1);
			}
		}
		genom.push_back(genomRow);
		genomRow.clear();

	}
	return genom;
}/* end of*/
 //同じ計算を何回もさせないためにNK地形を保存しておく。
vector<double_t> createEval()
{
	vector<double_t> evalList;
	evalList = comb_and_values(AllBitArrays(), calc_fit(NK_land(), IntMatrix(), AllBitArrays()));
	return evalList;
}

double_t maxEval(vector<double_t> evalList)
{
	double max;
	sort(evalList.begin(), evalList.end(), greater<double>());//降順ソート
	max = evalList[0];
	//cout << "max = " << max << endl;
	return max;
}

vector<vector<size_t>> select(vector<vector<size_t>> genome, vector<double_t> eval)
{
	vector<double_t> sortEval;
	vector<size_t> eliteRow;
	vector<vector<size_t>> elitegenome;

	sortEval = eval;
	sort(sortEval.begin(), sortEval.end(), std::greater<double>());//降順ソート

	for (int i = 0; i < MAX_GENOM_LIST; i++)
	{
		int a = -1;
		for (int j = 0; j < MAX_GENOM_LIST; j++)
		{
			if (sortEval[i] == eval[j])
			{
				//違う配列でも同じ評価値の場合
				a = j;
				for (int k = 0; k < N; k++)
				{
					eliteRow.push_back(genome[j][k]);
				}
				elitegenome.push_back(eliteRow);
				eliteRow.clear();
				break;
			}
		}
		eval[a] = -1;
	}

	for (size_t i = SELECT_GENOM; i < MAX_GENOM_LIST; i++)
	{
		elitegenome[i].erase(elitegenome[i].begin(), elitegenome[i].end());
	}

	return elitegenome;
}

vector<vector<size_t>> Genetic_recombination(vector<vector<size_t>> elitegenom,double_t Gene_combination_rate)
{
	vector<vector<size_t>> nextGeneration;
	vector<size_t> nextRows;
	vector<vector<size_t>> temp;
	vector<size_t> tempRow;

	vector<vector<size_t>> clone;
	int rand1 = 0;
	int onePoint;
	int twoPoint;


	//同じ個体の組み替えを2回行う
	for (size_t i = 0; i < SELECT_GENOM; i++)
	{
		if ((double)rand() / RAND_MAX <= Gene_combination_rate)
		{
			for (;;)
			{
				rand1 = rand() % SELECT_GENOM;
				if (rand1 != i)
				{
					break;
				}
			}
			onePoint = rand() % (N - 1);
			//cout << "1回目の"<< i << "個目の配列は" << rand1 << "個目の配列と二点交叉を行う" << endl;
			twoPoint = rand() % (N - onePoint - 1) + onePoint + 1;
			//cout << onePoint << "番目と" << twoPoint << "番目" << endl;

			for (int k = 0; k < N; k++)
			{
				if (k <= onePoint)
				{
					nextRows.push_back(elitegenom[i][k]);
				}
				else if (k > onePoint&& k < twoPoint)
				{
					nextRows.push_back(elitegenom[rand1][k]);
				}
				else
				{
					nextRows.push_back(elitegenom[i][k]);
				}
			}
			nextGeneration.push_back(nextRows);
			nextRows.clear();
		}
		else
		{
			for (int j = 0; j < N; j++)
			{
				nextRows.push_back(elitegenom[i][j]);
			}
			nextGeneration.push_back(nextRows);
			nextRows.clear();
		}
	}

	for (size_t i = 0; i < SELECT_GENOM; i++)
	{
		if ((double)rand() / RAND_MAX <= Gene_combination_rate)
		{
			for (;;)
			{
				rand1 = rand() % SELECT_GENOM;
				if (rand1 != i)
				{
					break;
				}
			}
			onePoint = rand() % (N - 1);
			//cout << "2回目の" <<i << "個目の配列は" << rand1 << "個目の配列と二点交叉を行う" << endl;
			twoPoint = rand() % (N - onePoint - 1) + onePoint + 1;
			//cout << onePoint << "番目と" << twoPoint << "番目" << endl;

			for (int k = 0; k < N; k++)
			{
				if (k <= onePoint)
				{
					nextRows.push_back(elitegenom[i][k]);
				}
				else if (k > onePoint&& k < twoPoint)
				{
					nextRows.push_back(elitegenom[rand1][k]);
				}
				else
				{
					nextRows.push_back(elitegenom[i][k]);
				}
			}
			nextGeneration.push_back(nextRows);
			nextRows.clear();
		}
		else
		{
			for (int j = 0; j < N; j++)
			{
				nextRows.push_back(elitegenom[i][j]);
			}
			nextGeneration.push_back(nextRows);
			nextRows.clear();
		}
	}

	return nextGeneration;
}
vector<vector<size_t>> mutation(vector<vector<size_t>> genome)
{
	vector<vector<size_t>> nextGenom;
	vector<size_t> nextRows;
	for (int i = 0; i < MAX_GENOM_LIST; i++) {
		if (INDIVIDUAL_MUTATION > (rand() % 100 + 1) / 100.0) {
			for (int j = 0; j < N; j++) {
				if (GENOM_MUTATION > (rand() % 100 + 1) / 100.0) {
					//cout << "突然変異するのは" << i << "番目のゲノム配列の" << j << endl;
					genome[i][j] == 0 ? (genome[i][j] = 1) : (genome[i][j] = 0);
				}
			}
		}
	}
	for (int i = 0; i < MAX_GENOM_LIST; i++)
	{
		for (int j = 0; j < N; j++)
		{
			nextRows.push_back(genome[i][j]);
		}
		nextGenom.push_back(nextRows);
		nextRows.clear();
	}

	return nextGenom;
}

vector<double_t> evaluation(vector<vector<size_t> > genome, vector<double_t> evalList)
{
	int total = 0;
	vector<double> x;
	for (int i = 0; i < MAX_GENOM_LIST; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (genome[i][j] == 1)
			{
				total += pow(2, j);
			}
		}
		x.push_back(evalList[total]);
		//cout << i << "番目の評価値は "<< x[i] << endl;
		total = 0;
	}
	return x;
}

vector<double_t> evaluation_MoveAvg(vector<vector<size_t>> genome, vector<double_t> eval, double_t translation_error_rate)
{
	size_t bitTotal = 0;
	double_t evalTotal = 0;
	vector<double_t> avg;
	vector<vector<size_t>> copy_genome;
	vector<size_t> copy_Rows;
	//0の時は移動平均なし
	if (Gene_product == 0)
	{
		for (size_t i = 0; i < MAX_GENOM_LIST; i++)
		{
			for (size_t j = 0; j < N; j++)
			{
				if (genome[i][j] == 1)
				{
					bitTotal += pow(2, j);
				}
			}
			evalTotal = eval[bitTotal];
			avg.push_back(evalTotal);
			bitTotal = 0;
			evalTotal = 0;
		}
	}
	else
	{
		for (size_t n = 0; n < MAX_GENOM_LIST; n++)
		{
			for (size_t i = 0; i < Gene_product; i++)
			{
				for (size_t j = 0; j < N; j++)
				{
					copy_Rows.push_back(genome[n][j]);
					if (translation_error_rate > (rand() % 100 + 1) / 100.0)
					{
						copy_Rows[j] == 0 ? (copy_Rows[j] = 1) : (copy_Rows[j] = 0);
					}

					if (copy_Rows[j] == 1)
					{
						bitTotal += pow(2, j);
					}

				}
				evalTotal += eval[bitTotal];
				bitTotal = 0;
				copy_genome.push_back(copy_Rows);
				copy_Rows.clear();
			}
			avg.push_back(evalTotal / Gene_product);
			copy_genome.clear();
			evalTotal = 0;
		}
	}

	return avg;
}

void printResult(vector<vector<size_t>> genom, int gn)
{
	cout << ": " << gn << "世代目 : " << endl;
	for (int i = 0; i < MAX_GENOM_LIST; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << genom[i][j];
		}
		cout << endl;
	}

}