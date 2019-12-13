/*
Program for finding frequencies of barcodes in FASTQ files of Chromium reads.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License.

Argument1: path to FASTQ file of reads #1 of Chromium paired end reads.
Output (STDOUT): frequencies of barcodes.

Compiling (requires g++ version >= 4.3, has been tested with gcc-8.2.0 on Ubuntu Linux):
g++ -std=c++11 chromium_reads_barcode_frequencies.cpp -o chromium_reads_barcode_frequencies
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <list>
#include <map>
using namespace std;


int main(int argc,  char **argv) {

	typedef map<string, int> string_int_map;
	string_int_map barcodes_map;

	if (argc == 2) {
		std::string file_path(argv[1]);
		std::ifstream file(file_path.c_str());
		if (file.good()) {
			std::string str;
			int counter = 0;
			int barcode_length = 16;

			while (std::getline(file, str)) {
				int size = str.length();
				int line_type = counter % 4;
				if (line_type == 1) {
					string barcode;
					if (size >= barcode_length) {
						barcode = str.substr(0, barcode_length);
						int barcode_count = barcodes_map.count(barcode) + 1;
						barcodes_map[barcode] = barcode_count;
					} else {
						cout << "Encountered a sequence that is too short" << endl;
						break;
						exit(EXIT_FAILURE);
					}
				}
				++counter;
			}
		} else {
			cout << file_path << ": file not found" << endl;
			exit(EXIT_FAILURE);
		}
		string_int_map::iterator pos;
		for (pos = barcodes_map.begin(); pos != barcodes_map.end(); ++pos) {
			cout << pos->first << "," << pos->second << endl;
		}
	} else {
		cout << argv[0] << ": tool for finding frequencies of barcodes in Chromium read FASTQ files" << endl;
		cout << "Usage: ";
		cout << argv[0] << " [FASTQ file]" << endl;
		exit(EXIT_FAILURE);
	}

    return 0;
}
