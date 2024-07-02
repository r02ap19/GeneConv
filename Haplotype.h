#pragma once

#include <map>
#include <set>
#include <stdio.h>
#include <iostream>
#include <fstream>

//maps containing whole sequence and markers: map<position, base> 
typedef std::map<int, char> bases;
typedef std::map<int, char> markers;

struct tract {
	int start; //Beginning of GC tract
	int length; //tract length (bp)
};

//map with true information on GC events: map<unique ID, tract>
typedef std::map<int, tract> conversions;

class Haplotype {
public:
	Haplotype();
	~Haplotype();

	bases sequence;
	markers SNPs;

	conversions GC_events;

	int number; // haplotype number

private:
};

