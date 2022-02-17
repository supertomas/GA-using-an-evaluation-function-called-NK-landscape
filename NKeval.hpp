#pragma once
#include<iostream>
#include<vector>
#include <time.h>
#include<random>
#include"math.h"
#include <chrono>

#define N 16
#define K 8
#define Landscape_seed 1


using namespace::std;
vector<vector<double>> NK_land();
vector<vector<size_t>> IntMatrix();
vector<vector<size_t>> AllBitArrays();
vector<vector<double_t>> calc_fit(vector<vector<double_t>>, vector<vector<size_t>>, vector<vector<size_t>>);
vector<double_t> comb_and_values(vector<vector<size_t>>, vector<vector<double_t>>);
