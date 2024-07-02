#pragma once

#define GenomeDK 0

#include <stdio.h>
#include <stdlib.h>
#if GenomeDK 
#include <unistd.h>
#else
#include <tchar.h> 
#include <direct.h>
#include <io.h>
#endif
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>

#include "Parameters.h"
#include "Haplotype.h"

using namespace std;

class Data {
public:
	Data();
	~Data();
	int holder;
	vector<Haplotype> haps;
	vector<Haplotype> tmp_haps;
	vector<Haplotype> converted_haps;
	vector<int> polymorphic_positions;

private:
};

clock_t extime;

Parameters para;
Data dataset;
string dir, dirOut;
ofstream seq, conv, summary, readData, joint;

std::random_device rd;
std::mt19937 rdgen(rd());

std::geometric_distribution<> geo(para.succes_p);
std::normal_distribution<> normal(para.mean_normal, para.var_normal);
std::uniform_int_distribution<> uni(para.min, para.max);
//std::_Beta_distribution<> beta(para.alpha, para.beta);
std::uniform_int_distribution<> WatsonCrick_bases(0, 3);

//Data representation
std::normal_distribution<> read_dist(para.mean_read_length, para.read_length_SD);


// Functions declaration
const string Int2Str(const int x);
void Simulate(void);
void initialize(void);
void add_SNPs(void);
void add_GeneConversions(int rep);
void outSeq_header(void);
void outConv_header(void);
void outSummary_header(void);
void outReadData_header(void);
void outJoint_header(void);
void outSeq(int h_number, char base, std::ofstream* out, int rep);
void outConv(int pos, int tract_length, int ID, std::ofstream* out, int rep);
vector<int> readEmpSNPs(const string& filename);
void outSummary(int rep, int genome_size, double detect_prob, std::ofstream* out, double class_1, double class_2, double class_3, double class_4, double class_5, double class_6, double class_7, double class_8, double class_9, double class_10, double class_n, double CO_detect);
void outReadData(int rep, int genome_size, double detect_prob, std::ofstream* out, double class_1, double class_2, double class_3, double class_4, double class_5, double class_6, double class_7, double class_8, double class_9, double class_10, double class_n, double CO_detect,
	double inconclusive_GC, double inconclusive_CO, double prob_nill, double prob_no_SNP_read);
void outJoint(int rep, int overlap_status, std::ofstream* out, int min_distance);
void post_pros(Data &dataset, int rep);

