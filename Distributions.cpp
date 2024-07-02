#include "Distributions.h"

#if GenomeDK
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path

	//para.SimNr = std::atoi(argv[1]);
	//para.succes_p = std::atof(argv[2]);

	Simulate();
	cout << "Simulation completed" << endl;
	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[])
{
	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path
	extime = clock();

	Simulate();
	std::cout << "Simulation completed" << endl;
	return 0;
}
#endif
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}

void Simulate(void) {
	std::cout << "Simulation nr. " << para.SimNr << endl;
	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Para.txt";

	outSeq_header();
	outConv_header();
	outSummary_header();
	outReadData_header();
	outJoint_header();

	for (int r = 0; r < para.reps; r++) {

		cout << "Replicate: " << r << " out of " << para.reps << "\n" << endl;

		cout << "Initializing...";
		initialize();
		cout << "done" << endl;

		extime = clock() - extime;
		std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
		extime = clock();
		cout << "\n" << endl;

		cout << "Adding markers...";
		add_SNPs();
		cout << "done" << endl;

		extime = clock() - extime;
		std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
		extime = clock();
		cout << "\n" << endl;

		cout << "Conducting gene conversion events...";
		add_GeneConversions(r);
		cout << "done" << endl;

		extime = clock() - extime;
		std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
		extime = clock();
		cout << "\n" << endl;

		if (para.post_processing == true) {
			cout << "Conducting post-processing..." << endl;
			post_pros(dataset, r);
			cout << "done" << endl;

			extime = clock() - extime;
			std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
			extime = clock();
			cout << "\n" << endl;
		}
	}

	para.outPara(name);
	if (seq.is_open()) seq.close();
	if (conv.is_open()) conv.close();
	if (summary.is_open()) summary.close();
	if (readData.is_open()) readData.close();
	if (joint.is_open()) joint.close();

}

void initialize(void) {
	Haplotype H1, H2;
	int base_value, base_counter, distance;
	std::map<int, char>::iterator iter;
	base_counter = 0;
	std::gamma_distribution<> dists(para.shape, para.rate);

	if (para.coalescent_dist == true) {
		int N_SNPs = round(para.SNP_density * para.bp);
		for (int x = 0; x < N_SNPs; x++) {
			distance = round(dists(rdgen));
			if (distance < 1) {
				distance = 1;
			}
			base_counter = base_counter + distance;
			dataset.polymorphic_positions.push_back(base_counter);
		}
		para.realized_genome_size = base_counter;
	} else {
		para.realized_genome_size = para.bp;
	}

	if (para.empirical_dist == true) {
		string filename = para.data_name;
		vector<int> distances = readEmpSNPs(filename);
		int sum = 0;
		for (int x = 0; x < distances.size(); x++) {
			sum += distances[x];
		}
		para.realized_genome_size = sum;
		para.SNP_density = double(distances.size()) / double(sum);
		distances.clear();
	}

	for (int x = 0; x < para.realized_genome_size; x++) {
		base_value = WatsonCrick_bases(rdgen);
		if (base_value == 0)(H1.sequence[x] = 'A');
		if (base_value == 1)(H1.sequence[x] = 'T');
		if (base_value == 2)(H1.sequence[x] = 'C');
		if (base_value == 3)(H1.sequence[x] = 'G');
	}

	H2 = H1;
	H2.number = 2;
	dataset.tmp_haps.push_back(H1);
	dataset.tmp_haps.push_back(H2);
	H1.sequence.clear();
	H2.sequence.clear();
}

void add_SNPs(void) {
	int N_SNPs = round(para.SNP_density * para.bp);
	int pos;
	double prob;
	std::map<int, char>::iterator iter;
	std::vector<Haplotype>::iterator iter2;

	std::uniform_int_distribution<> genome_pos(0, para.realized_genome_size);
	std::uniform_real_distribution<> SNP_prob(0.0, 1.0);

	if (para.coalescent_dist == true) {
		for (int x = 0; x < (N_SNPs-1); x++) {
			dataset.tmp_haps[0].sequence[dataset.polymorphic_positions[x]] = 'X';
			dataset.tmp_haps[0].SNPs[dataset.polymorphic_positions[x]] = 'X';
		}
	}
	if (para.binom_dist == true) {
		for (int x = 0; x < para.realized_genome_size; x++) {
			prob = SNP_prob(rdgen);
			if (prob < para.p_success) {
				dataset.tmp_haps[0].sequence[x] = 'X';
				dataset.tmp_haps[0].SNPs[x] = 'X';
			}
		}
	} 

	if (para.empirical_dist == true) {
		string filename = para.data_name;
		vector<int> distances = readEmpSNPs(filename);
		int cumulative_dist = 0;
		for (int x = 0; x < (distances.size()-1); x++) {
			cumulative_dist = cumulative_dist + distances[x];
			dataset.tmp_haps[0].sequence[cumulative_dist] = 'X';
			dataset.tmp_haps[0].SNPs[cumulative_dist] = 'X';
		}
		distances.clear();
	}

	if (para.binom_dist == false && para.coalescent_dist == false && para.empirical_dist == false) {
		pos = genome_pos(rdgen);
		for (int x = 0; x < N_SNPs; x++) {
			pos = genome_pos(rdgen);
			while (dataset.tmp_haps[0].sequence[pos] == 'X') {
				pos = genome_pos(rdgen);
			}
			dataset.tmp_haps[0].sequence[pos] = 'X';
			dataset.tmp_haps[0].SNPs[pos] = 'X';
		}
	}
	if (!dataset.polymorphic_positions.empty())(dataset.polymorphic_positions.clear());

	dataset.haps = dataset.tmp_haps;
	dataset.tmp_haps.clear();

}

//GC function
void add_GeneConversions(int rep) {
	int pos, GC_tractLength;
	double peneVal;
	bool sample;
	std::map<int, char>::iterator iter;
	std::map<int, tract>::iterator iter2;
	//std::map<int, tract>::iterator iter3;
	tract event;
	std::uniform_int_distribution<> genome_pos(0, para.realized_genome_size);
	std::geometric_distribution<> geo(para.succes_p);
	std::uniform_real_distribution<> pene(0.0, 1.0);
	std::vector<int> lengths;
	std::vector<int> GC_positions;

	sample = true;

	while (sample == true) {
		sample = false;
		lengths.clear();
		GC_positions.clear();
		for (int x = 0; x < para.GC_events; x++) {
			pos = genome_pos(rdgen);
			GC_positions.push_back(pos);
		}
		sort(GC_positions.begin(), GC_positions.end());

		for (int x = 0; x < para.GC_events; x++) {
			if (para.distribution == "geometric") {
				GC_tractLength = (geo(rdgen) + 1); //converting geometric dist to count number of trials before success, rather than number of failures before success (see C++ documentation). 
			}
			if (para.distribution == "normal") {
				GC_tractLength = round(normal(rdgen));
			}
			if (para.distribution == "uniform") {
				GC_tractLength = round(uni(rdgen));
			}
			lengths.push_back(GC_tractLength);
		}
		for (int x = 0; x < para.GC_events; x++) {
			if ((GC_positions[x] + lengths[x]) > para.realized_genome_size) {
				sample = true;
			}
		}
	}
	for (int x = 0; x < para.GC_events; x++) {
		event.start = GC_positions[x];
		event.length = lengths[x];
		dataset.haps[1].GC_events[x] = event;
		peneVal = 1.0;
		for (int y = 0; y < lengths[x]; y++) {
			dataset.haps[1].sequence[GC_positions[x] + y] = dataset.haps[0].sequence[GC_positions[x] + y];
			if (dataset.haps[0].sequence[GC_positions[x] + y] == 'X') {
				peneVal = pene(rdgen);
				if (peneVal < para.penetrance) {
					dataset.haps[1].SNPs[GC_positions[x] + y] = 'X';
				}
			}
		}
	}

	lengths.clear();
	GC_positions.clear();

	if (para.post_processing == false) {
		cout << "Writing output files...";

		dataset.converted_haps = dataset.haps;

		for (iter2 = dataset.converted_haps[1].GC_events.begin(); iter2 != dataset.converted_haps[1].GC_events.end(); iter2++) {
			outConv(iter2->second.start, iter2->second.length, iter2->first, &conv, rep);
		}
		for (iter = dataset.converted_haps[0].sequence.begin(); iter != dataset.converted_haps[0].sequence.end(); iter++) {
			outSeq(dataset.converted_haps[0].number, iter->second, &seq, rep);
		}
		for (iter = dataset.converted_haps[1].sequence.begin(); iter != dataset.converted_haps[1].sequence.end(); iter++) {
			outSeq(dataset.converted_haps[1].number, iter->second, &seq, rep);
		}

		dataset.converted_haps.clear();
		dataset.haps.clear();

		cout << "done" << endl;

		extime = clock() - extime;
		std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
		extime = clock();
		cout << "\n" << endl;

	}
}

Data::Data() {
	holder = 0;
}
Data::~Data() {
}

vector<int> readEmpSNPs(const string& filename) {
	ifstream infile(filename);
	cout << "Name of the data input file: " << filename << endl;

	if (!infile.is_open()) {
		cerr << "Error opening file containing empirical inter SNP distances; file name: " << filename << std::endl;
		return std::vector<int>();
	}

	vector<int> distances;
	int distance;

	while (infile >> distance) {
		distances.push_back(distance);
	}

	infile.close();
	return distances;
}

void post_pros(Data &dataset, int rep) {

	map<int, char>::iterator iter;
	int counter = 0, cross_pos = 0, read_length = 0, read_start = 0, read_end = 0, SNP_count_forward = 0, SNP_count_backward = 0, SNP_count = 0;
	double counter2, detect_prob, detect_prob_CO = 0.0, CO_count = 0.0, inconclusive_event_CO = 0.0, inconclusive_event_GC = 0.0,
		no_SNP_read = 0.0, nill = 0.0, prob_no_SNP_read = 0.0, prob_nill= 0.0, read_scaler = 0.0;
	int pos_plus = 0, count_1 = 0, count_2 = 0, SNP_pos_1 = 0, SNP_pos_2 = 0;
	vector<double> proportions;
	vector<int> conversion_events;
	vector<int> cross_positions;
	vector<int> overlap_status;
	vector<int> joint_tract_distances;

	//GC calling, full reads
	for (iter = dataset.haps[0].SNPs.begin(); iter != dataset.haps[0].SNPs.end(); iter++) {
		pos_plus = 0;
		while (dataset.haps[1].SNPs[iter->first + pos_plus] == 'X') {
			counter += 1;
			pos_plus = ((next(iter, counter)->first) - (iter->first));
		}
		if (pos_plus < 0 || pos_plus > para.realized_genome_size) {
			break;
		}
		if (counter > 0) {
			conversion_events.push_back(counter);
		}
		if (!(next(iter, counter) == dataset.haps[0].SNPs.end())){
			advance(iter, counter);
		}
		counter = 0;
	}

	counter2 = 0.0;
	for (int x = 1; x < 11; x++) {
		proportions.push_back((count(conversion_events.begin(), conversion_events.end(), x) / fmax(conversion_events.size(),1.0)));
		counter2 += count(conversion_events.begin(), conversion_events.end(), x);
	}
	if (!conversion_events.empty()) {
		if ((*max_element(conversion_events.begin(), conversion_events.end()) > 10)) {
			proportions.push_back((conversion_events.size() - counter2) / conversion_events.size());
		}
	}
	detect_prob = (conversion_events.size() / (double)para.GC_events);

	//Calculate detection of CO events
	//Sample CO positions    
	CO_count = 0;
	inconclusive_event_CO = 0.0;
	inconclusive_event_GC = 0.0;
	for (int i = 0; i < para.n_cross; i++) {
		//place and build read
		std::uniform_int_distribution<> genome_pos(0, para.realized_genome_size);
		read_start = genome_pos(rdgen);
		read_end = read_start + round(read_dist(rdgen));

		if (read_end > para.realized_genome_size)(read_end = para.realized_genome_size);

		std::uniform_real_distribution<> scaled_reads(0, 1); //use to scale prob of sampling reads (the longer the read, the more likely to contain a crossover event
		read_scaler = scaled_reads(rdgen);

		while ((((read_end - read_start) - para.min_read) < 1) || ((para.max_read - (read_end - read_start)) < 1) || ((double(read_end - read_start) / double(para.max_read)) < read_scaler)) { //check read isn't too long, too short, and scale prob of containing crossover by read length
			read_start = genome_pos(rdgen);
			read_end = read_start + round(read_dist(rdgen));
			read_scaler = scaled_reads(rdgen);
		}


		//sample CO event
		std::uniform_int_distribution<> crossover_position(read_start, read_end);
		cross_pos = crossover_position(rdgen);
		//Check if crossover is callable 
		SNP_count_forward = 0;
		SNP_count_backward = 0;
		for (int j = cross_pos; j < read_end + 1; j++) {
			if (dataset.haps[0].sequence[j] == 'X') {
				SNP_count_forward += 1;
			}
		}
		for (int j = read_start; j < cross_pos; j++) {
			if (dataset.haps[0].sequence[j] == 'X') {
				SNP_count_backward += 1;
			}
		}
		if (SNP_count_forward > 1 && SNP_count_backward > 1) {
			CO_count += 1;
		}
		if ((SNP_count_forward == 1 && SNP_count_backward > 1) || (SNP_count_forward > 1 && SNP_count_backward == 1) || (SNP_count_forward == 1 && SNP_count_backward == 1)) {
			inconclusive_event_CO += 1.0;
		}
		if (SNP_count_forward == 0 || SNP_count_backward == 0) {
			nill += 1.0;
		}
		if (SNP_count_forward == 0 && SNP_count_backward == 0){
			no_SNP_read += 1.0;
		}

	}
	detect_prob_CO = (double)CO_count / (double)para.n_cross;
	prob_nill = nill/(double)para.n_cross;
	prob_no_SNP_read = no_SNP_read / (double)para.n_cross;

	if (proportions.size() > 10) {
		outSummary(rep, para.realized_genome_size, detect_prob, &summary, proportions[0], proportions[1], proportions[2], proportions[3], proportions[4], proportions[5], proportions[6], proportions[7], proportions[8], proportions[9], proportions[10], detect_prob_CO);
	}
	else {
		outSummary(rep, para.realized_genome_size, detect_prob, &summary, proportions[0], proportions[1], proportions[2], proportions[3], proportions[4], proportions[5], proportions[6], proportions[7], proportions[8], proportions[9], 0.0, detect_prob_CO);
	}


	conversion_events.clear();
	proportions.clear();

	detect_prob = 0.0;

	//GC calling, sampled read length
	read_start = 0;
	read_end = read_start + round(read_dist(rdgen));
	for (iter = dataset.haps[0].SNPs.begin(); iter != dataset.haps[0].SNPs.end(); iter++) {

		if (read_end > para.realized_genome_size)(read_end = para.realized_genome_size);

		while ((iter->first > read_end) || (((read_end - read_start) - para.min_read) < 1)){
			read_start = read_end;
			read_end += round(read_dist(rdgen));
		}

		pos_plus = 0;
		while (dataset.haps[1].SNPs[iter->first + pos_plus] == 'X') {
			counter += 1;
			pos_plus = ((next(iter, counter)->first) - (iter->first));
		}
		if (pos_plus < 0 || pos_plus > para.realized_genome_size) {
			break;
		}
		if (counter > 0) {
			SNP_count_forward = 0;
			SNP_count_backward = 0;
			for (int j = iter->first + pos_plus + 1; j < read_end + 1; j++) {
				if (dataset.haps[0].sequence[j] == 'X' && dataset.haps[1].sequence[j] != 'X') {
					SNP_count_forward += 1;
				}
				if (SNP_count_forward > 1) {
					break;
				}
			}

			for (int j = read_start; j < iter->first + pos_plus + 1; j++) {
				if (dataset.haps[0].sequence[j] == 'X' && dataset.haps[1].sequence[j] != 'X') {
					SNP_count_backward += 1;
				}
				if (SNP_count_backward > 1) {
					break;
				}
			}
			if ((SNP_count_forward == 0 && SNP_count_backward > 0) || (SNP_count_forward > 0 && SNP_count_backward == 0)) {
				inconclusive_event_GC += 1.0;
			}

			if (SNP_count_forward > 0 && SNP_count_backward > 0) {
				conversion_events.push_back(counter);
			}
		}
		if (!(next(iter, counter) == dataset.haps[0].SNPs.end())) {
			advance(iter, counter);
		}
		counter = 0;
	}

	counter2 = 0.0;
	for (int x = 1; x < 11; x++) {
		proportions.push_back((count(conversion_events.begin(), conversion_events.end(), x) / fmax(conversion_events.size(), 1.0)));
		counter2 += count(conversion_events.begin(), conversion_events.end(), x);
	}
	if (!conversion_events.empty()) {
		if ((*max_element(conversion_events.begin(), conversion_events.end()) > 10)) {
			proportions.push_back((conversion_events.size() - counter2) / conversion_events.size());
		}
	}
	detect_prob = (conversion_events.size() / (double)para.GC_events);
	inconclusive_event_CO = (inconclusive_event_CO / ((double)para.n_cross));
	inconclusive_event_GC = (inconclusive_event_GC / ((double)para.GC_events));

	if (proportions.size() > 10) {
		outReadData(rep, para.realized_genome_size, detect_prob, &readData, proportions[0], proportions[1], proportions[2], proportions[3], proportions[4], proportions[5], proportions[6], proportions[7], proportions[8], proportions[9], proportions[10], detect_prob_CO, inconclusive_event_GC, inconclusive_event_CO, prob_nill, prob_no_SNP_read);
	}
	else {
		outReadData(rep, para.realized_genome_size, detect_prob, &readData, proportions[0], proportions[1], proportions[2], proportions[3], proportions[4], proportions[5], proportions[6], proportions[7], proportions[8], proportions[9], 0.0, detect_prob_CO, inconclusive_event_GC, inconclusive_event_CO, prob_nill, prob_no_SNP_read);
	}
	if (para.joint_tracts == true) {
		for (int i = 0; i < dataset.haps[1].GC_events.size() - 1; ++i) {
			count_1 = 0;
			SNP_pos_1 = 0;

			for (int j = 0; j < dataset.haps[1].GC_events[i].length; ++j) {
				if (dataset.haps[1].sequence[(dataset.haps[1].GC_events[i].start + j)] == 'X') {
					count_1 = count_1 + 1;
					if (SNP_pos_1 == 0)(SNP_pos_1 = (dataset.haps[1].GC_events[i].start + j));
				}
			}
			count_2 = 0;
			SNP_pos_2 = 0;
			for (int j = 0; j < dataset.haps[1].GC_events[i+1].length; ++j) {
				if (dataset.haps[1].sequence[(dataset.haps[1].GC_events[i+1].start + j)] == 'X') {
					count_2 = count_2 + 1;
					SNP_pos_2 = (dataset.haps[1].GC_events[i+1].start + j);
				}
			}

			int SNP_count = 0;
			if ((count_1 > 0 && count_2 > 0)) {
				for (int j = 0; j < (((dataset.haps[1].GC_events[i+1].start + dataset.haps[1].GC_events[i + 1].length) - dataset.haps[1].GC_events[i].start)); ++j) {
					if (dataset.haps[0].sequence[((dataset.haps[0].GC_events[i].start) + j)] == 'X') {
						SNP_count = SNP_count + 1;
					}
				}
				if ((count_1 + count_2) == SNP_count) {
					overlap_status.push_back(1);
					joint_tract_distances.push_back(SNP_pos_2 - SNP_pos_1);
				} else {
					overlap_status.push_back(0);
					joint_tract_distances.push_back(-99);
				}
			} else {
				overlap_status.push_back(0);
				joint_tract_distances.push_back(-99);
			}
		}
		for (int x = 0; x < dataset.haps[1].GC_events.size() - 1; x++) {
			outJoint(rep, overlap_status[x], &joint, joint_tract_distances[x]);
		}
	}

	conversion_events.clear();
	proportions.clear();
	overlap_status.clear();
	joint_tract_distances.clear();
	dataset.haps.clear();
}

void outSeq_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outSeq.txt";
	seq.open(name.c_str());
	seq << "Haplotype\tseq\trep" << endl;
}

void outConv_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outConv.txt";
	conv.open(name.c_str());
	conv << "pos\ttract_length\tID\trep" << endl;
}

void outSummary_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outSummary.txt";
	summary.open(name.c_str());
	summary << "rep\tgenome_size\tdetect_prob\tclass_1\tclass_2\tclass_3\tclass_4\tclass_5\tclass_6\tclass_7\tclass_8\tclass_9\tclass_10\tclass_n\tCO_detect" << endl;
}

void outReadData_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outReadData.txt";
	readData.open(name.c_str());
	readData << "rep\tgenome_size\tdetect_prob\tclass_1\tclass_2\tclass_3\tclass_4\tclass_5\tclass_6\tclass_7\tclass_8\tclass_9\tclass_10\tclass_n\tCO_detect\tIF_GC\tIF_CO\tprob_nill\tprob_no_SNP_read" << endl;
}

void outJoint_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outJoint.txt";
	joint.open(name.c_str());
	joint << "rep\toverlap_status\tmin_distance\t" << endl;
}

void outSeq(int h_number, char base, std::ofstream* out, int rep) {
	*out << h_number << "\t" << base << "\t" << rep << "\t";
	*out << endl;
}

void outConv(int pos, int tract_length, int ID, std::ofstream* out, int rep) {
	*out << pos << "\t" << tract_length << "\t" << ID << "\t" << rep << "\t";
	*out << endl;
}

void outSummary(int rep, int genome_size, double detect_prob, std::ofstream* out, double class_1, double class_2, double class_3, double class_4, double class_5, double class_6, double class_7, double class_8, double class_9, double class_10, double class_n, double CO_detect) {
	*out << rep << "\t" << genome_size << "\t" << detect_prob << "\t" << class_1 << "\t" << class_2 << "\t" << class_3 << "\t" << class_4 << "\t" << class_5 << "\t" << class_6 << "\t" << class_7 << "\t" << class_8 << "\t" << class_9 << "\t" << class_10 << "\t" << class_n << "\t" << CO_detect << "\t";
	*out << endl;
}

void outReadData(int rep, int genome_size, double detect_prob, std::ofstream* out, double class_1, double class_2, double class_3, double class_4, double class_5, double class_6, double class_7, double class_8, double class_9, double class_10, double class_n, double CO_detect, double inconclusive_GC, double inconclusive_CO, double prob_nill, double prob_no_SNP_read){
	*out << rep << "\t" << genome_size << "\t" << detect_prob << "\t" << class_1 << "\t" << class_2 << "\t" << class_3 << "\t" << class_4 << "\t" << class_5 << "\t" << class_6 << "\t" << class_7 << "\t" << class_8 << "\t" << class_9 << "\t" << class_10 << "\t" << class_n << "\t" << CO_detect << "\t" << inconclusive_GC << "\t" << inconclusive_CO << "\t" << prob_nill << "\t" << prob_no_SNP_read << "\t";
	*out << endl;
}

void outJoint(int rep, int overlap_status, std::ofstream* out, int min_distance) {
	*out << rep << "\t" << overlap_status << "\t" << min_distance << "\t";
	*out << endl;
}