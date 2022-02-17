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
	globalEvalList = createEval();//���̏����Ɏ��Ԃ�������B
	global_max_eval = maxEval(globalEvalList);
	ofs << ",,,,,�|��G���[��" << endl;
	for (size_t i = 0; i < 11; i++)
	{
		ofs << "," << i * 0.02;
	}
	ofs << endl;
	cout << "���s�񐔂���͂��Ă��������B" << endl;
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
						movegenome.genomeArray = Genetic_recombination(select(movegenome.genomeArray, movegenome.eval), Gene_combination_rate);//��`�q��g�݊�����B
						movegenome.genomeArray = mutation(movegenome.genomeArray);//�ˑR�ψق�����B
					}
					else
					{

						movegenome.eval = evaluation(movegenome.genomeArray, globalEvalList);
						generation_Eval = movegenome.eval;//�������g�݊������邽�߂ɕʂ�vector��genom.eval��ۑ����Ă����B
						sort(generation_Eval.begin(), generation_Eval.end(), greater<double_t>());//�~���\�[�g
						total += generation_Eval[0] / global_max_eval;
						cout << "Recom_rate : " << Gene_combination_rate << endl;
						cout << "translation error" << translation_error_rate << endl;
						cout << i << "��� : " << generation_Eval[0] / global_max_eval << endl;
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
	ofs << "�n�`�̗�����" << Landscape_seed << endl;
	ofs << "N," << N << ",K," << K << endl;
	ofs << "�S�̂̌̐�," << MAX_GENOM_LIST << ",�G���[�g�̐�," << SELECT_GENOM << ",���㐔," << MAX_GENERATION << endl;
	ofs << "��`�q�Y����," << Gene_product << ",��`�q�ψٗ�," << GENOM_MUTATION << endl;
}