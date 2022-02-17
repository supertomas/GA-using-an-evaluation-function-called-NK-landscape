#include"NKeval.hpp"
#define MAX_GENOM_LIST 20
#define SELECT_GENOM 10
#define INDIVIDUAL_MUTATION 1.0
#define GENOM_MUTATION 0.04
#define MAX_GENERATION 20
#define step_mutation 0.02
#define step_recom 0.05
#define Gene_product 80

struct moveGenome
{
	vector<vector<size_t>> genomeArray;
	vector<double_t> eval;
};

vector<vector<size_t>> init();
double_t maxEval(vector<double_t>);
vector<double_t> createEval();
vector<double_t> evaluation(vector<vector<size_t> >, vector<double_t>);
vector<double_t> evaluation_MoveAvg(vector<vector<size_t> >, vector<double_t>, double_t);
vector<size_t> Conversation_Decimal(vector<vector<size_t>>);
vector<vector<size_t>> select(vector<vector<size_t>>, vector<double_t>);
vector<vector<size_t>> Genetic_recombination(vector<vector<size_t>>,double_t);
vector<vector<size_t>> mutation(vector<vector<size_t>>);
void printResult(vector<vector<size_t>>, int);