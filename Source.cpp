#pragma once
#include"GAalgorithm.hpp"
#include <fstream>
#include "time.h"
#define genomeseed 20
int main()
{
	moveGenome movegenome;

	vector<double_t> globalEvalList, generation_Eval;

	double_t global_max_eval;
	time_t now = time(nullptr);
	const char* fileName = "MaxevalList.csv";
	ofstream ofs(fileName);
	double_t translation_error_rate = 0.0, eval_avg = 0.0, total = 0.0, Gene_combination_rate = 0.05;
	size_t trials;
	vector<vector<size_t>> Save_first_genomes;
	vector<vector<vector<size_t>>> Save_first_genomes_list;

	srand(Landscape_seed);
	globalEvalList = createEval();//この処理に時間がかかる。
	global_max_eval = maxEval(globalEvalList);
	ofs << ",,,,,翻訳エラー率" << endl;
	for (size_t i = 0; i < 11; i++)
	{
		ofs << "," << i * 0.02;
	}
	ofs << endl;
	cout << "試行回数を入力してください。" << endl;
	cin >> trials;
	srand(genomeseed);
	for (size_t i = 0; i < trials; i++)
	{
		Save_first_genomes = init();
		Save_first_genomes_list.push_back(Save_first_genomes);
		Save_first_genomes.clear();
	}

	while (Gene_combination_rate <= 0.45)
	{
		ofs << Gene_combination_rate;
		while (translation_error_rate <= 0.2)
		{
			for (size_t i = 0; i < trials; i++)
			{
				movegenome.genomeArray = Save_first_genomes_list[i];
				for (size_t n = 1; n <= MAX_GENERATION; n++)
				{
					if (n != MAX_GENERATION)
					{
						movegenome.eval = evaluation_MoveAvg(movegenome.genomeArray, globalEvalList, translation_error_rate);
						movegenome.genomeArray = Genetic_recombination(select(movegenome.genomeArray, movegenome.eval), Gene_combination_rate);//遺伝子を組み換える。
						movegenome.genomeArray = mutation(movegenome.genomeArray);//突然変異させる。
					}
					else
					{

						movegenome.eval = evaluation(movegenome.genomeArray, globalEvalList);
						generation_Eval = movegenome.eval;//正しく組み換えするために別のvectorにgenom.evalを保存しておく。
						sort(generation_Eval.begin(), generation_Eval.end(), greater<double_t>());//降順ソート
						total += generation_Eval[0] / global_max_eval;
						cout << "Recom_rate : " << Gene_combination_rate << endl;
						cout << "translation error" << translation_error_rate << endl;
						cout << i << "回目 : " << generation_Eval[0] / global_max_eval << endl;
						generation_Eval.clear();
					}
				}
			}
			eval_avg = total / trials;
			ofs << "," << eval_avg;
			movegenome.genomeArray.clear();
			translation_error_rate += step_mutation;
			total = 0.0;
		}
		Gene_combination_rate += step_recom;
		translation_error_rate = 0.0;
		ofs << endl;
	}
	ofs << "地形の乱数種" << Landscape_seed << endl;
	ofs << "N," << N << ",K," << K << endl;
	ofs << "全体の個体数," << MAX_GENOM_LIST << ",エリート個体数," << SELECT_GENOM << ",世代数," << MAX_GENERATION << endl;
	ofs << "遺伝子産物数," << Gene_product << ",遺伝子変異率," << GENOM_MUTATION << endl;
}