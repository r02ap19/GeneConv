#pragma once

#include <fstream>

using namespace std;

class Parameters
{
public:
	Parameters();
	~Parameters();
	//Simulation
	int SimNr;
	bool post_processing;
	int reps; // Number of simulation replicates to run
	bool joint_tracts; //if true, the occurence of multiple tracts appearing in data as a single tract is checked.

	//Genome
	int bp; //Number of base pairs in haplotype
	double SNP_density; //Density of SNPs
	double rate; //rate parameter for gamma distribution of inter SNP distances (see Field et al. 2016).
	double shape; //shape parameter for gamma distribution of inter SNP distances (see Field et al. 2016).
	bool coalescent_dist; // if true, distances between SNPs are gamma distributed, approximating the expectation from coalescent theory (see Field et al. 2016).
	bool empirical_dist; //if true, a file with interSNP distances from HiFi PacBio data is loaded and used to distribute SNPs.
	int realized_genome_size; //size of the genome when genome is is determined by inter-SNP distances i.e. when coalescent_dist == true. Else, realized_genome_size = bp.
	bool binom_dist; //if true, each base has a specific probability of being a SNP (p_success) such that the distribution of SNPs (and inter-SNP distances) follow a binomial dist.
	double p_success; //per base probability of being a SNP
	string data_name; //name of the txt file with empirical data on inter SNP distances if empirical_cist == true. 

	//Gene conversions (GC)
	int GC_events; // number of gene conversion events (fixed for now, but could be sampled).
	string distribution; // distribution of GC locations: geometric, normal, uniform, beta or exponential
	double succes_p; // probability of success, geometric distribution
	double mean_normal; // mean of normal distribution
	double var_normal; // variance, normal distribution
	int min; //lower bound for uniformly distributed GC tract lengths
	int max; //upper bound for uniformaly distributed GC tract lengths
	double alpha; //alpha parameter in beta distribution for beta distributed GC tract lengths
	double beta; //beta parameter in beta distribution for beta distributed GC tract lengths
	double penetrance; //probability of converting a SNP given that it is in a NCO tract.

	//Crossovers (CO)
	int n_cross; // number of crossover events for the entire genome

	//Data representation
	int mean_read_length; //mean length of reads (sampled from normal distribution)
	int read_length_SD; //standard deviation in read length (sampled from normal distribution)
	int min_read; //lower bound for read size.
	int max_read; //upper bound for read size.

	void outPara(string dir); //Parameter output
};


