#include "Parameters.h"

Parameters::Parameters() {

	//Simulation
	SimNr = 1;
	post_processing = true;
	reps = 1000;
	joint_tracts = false;

	//Genome
	bp = 500000000;
	SNP_density = 0.00083; // given parameters below, acheived approximate SNP density of Zhao et al. 2003: 0.00083.  from old dataset: 0.00059375
	shape = 0.25567091;
	rate = 237.3453398;// scale parameter equivalent to the rate parameter 0.00421327 inferred by Peter from the PacBio data//OLD:scale parameter equivalent to the rate parameter 1/0.00166 yielding approx 0.00083 SNPs/base and 10Mbp genome
	coalescent_dist = false;
	empirical_dist = false;
	binom_dist = true;
	p_success = 0.00083;
	data_name = "1.txt";

	//Gene conversions (GC)
	GC_events = 5000;
	distribution = "geometric"; // distribution of GC lengths: geometric, normal, uniform beta or exponential
	succes_p = (1.0/50.0);
	mean_normal = 100.0;
	var_normal = 10;
	min = 1;
	max = 10;
	alpha = 0.5;
	beta = 0.5;
	penetrance = 1.0;

	//Crossovers (CO)
	n_cross = 1000;

	//Data representation
	mean_read_length = 12256;
	read_length_SD = 7484;
	min_read = 3314;
	max_read = 26356;
}

Parameters::~Parameters() {
}

void Parameters::outPara(string name) {

	ofstream out;
	out.open(name.c_str());

	//Simulation
	out << "SimNr\t" << SimNr << endl;

	//Genome
	out << "bp\t" << bp << endl;
	out << "SNP_density\t" << SNP_density << endl;
	out << "rate\t" << rate << endl;
	out << "shape\t" << shape << endl;
	out << "realized_genome_size\t" << realized_genome_size << endl;
	out << "coalescent_dist\t" << coalescent_dist << endl;
	out << "empirical_dist\t" << empirical_dist << endl;
	out << "binom_dist\t" << binom_dist << endl;
	out << "p_success\t" << p_success << endl;
	out << "data_name\t" << data_name << endl;

	//Gene conversions
	out << "GC_events\t" << GC_events << endl;
	out << "distribution\t" << distribution << endl;
	out << "succes_p\t" << succes_p << endl;
	out << "mean_normal\t" << mean_normal << endl;
	out << "var_normal\t" << var_normal << endl;
	out << "min\t" << min << endl;
	out << "max\t" << max << endl;
	out << "alpha\t" << alpha << endl;
	out << "beta\t" << beta << endl;

	out << "post_processing\t" << post_processing << endl;
	out << "replicates\t" << reps << endl;

	//Crossovers (CO)
	out << "n_cross\t" << n_cross << endl;

	//data representation
	out << "mean_read_length\t" << mean_read_length << endl;
	out << "read_length_SD\t" << read_length_SD << endl;
	out << "min_read\t" << min_read << endl;
	out << "joint_tracts\t" << joint_tracts << endl;
	out << "penetrance\t" << penetrance << endl;

	out.close();
}

