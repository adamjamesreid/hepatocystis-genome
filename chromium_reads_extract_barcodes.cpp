/*
C++ program for extracting Chromium barcode sequences of Chromium reads that have been selected for assembly.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License.

The aim of the pipeline that this program is a part of is to sort Chromium reads into assembly batches with a fixed maximum number of reads, where in each assembly batch as many reads as possible share the same barcode.
In order to do this, a table is made where the first column contains barcode sequences and the second column contains read names (this is the output of this program).

After sorting the output table alphabetically by the barcode sequences and then splitting the table into smaller tables with a fixed number of rows, 
each of the smaller tables will correspond to an assembly batch. Column 2 will contain the names of reads, so the selected reads can then be extracted from FASTQ files by their names files using SeqTK.

Argument1: path to a text file with read names
Argument2: path to reads FASTA file, derived from FASTQ file (#1 file of paired end reads)
Argument3: path to a CSV file, with the first column containing sequences of barcodes that have more than one read associated with them
	The header line will be skipped by this program
Output: (printed to STDOUT)
	first column: read barcode (extracted from the FASTA file), second column: read name.
	If the barcode is not listed in the selected barcodes file (argument3), "NA" is printed in place of barcode

Compiling (requires g++ version >= 4.3, has been tested with gcc-8.2.0 on Ubuntu Linux):
g++ -std=c++11 chromium_reads_extract_barcodes.cpp -o chromium_reads_extract_barcodes
 */

#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <cstdio>
#include <cerrno>
#include <sstream>
#include <stdlib.h>
#include <set>
using namespace std;


string get_file_contents(const char *filename) {
	FILE *fp = fopen(filename, "rb");
	if (fp) {
		string contents;
		fseek(fp, 0, SEEK_END);
		contents.resize(ftell(fp));
		rewind(fp);
		fread(&contents[0], 1, contents.size(), fp);
		fclose(fp);
		return(contents);
	}
	throw(errno);
}


vector<string> split(const string &s, char delim) {
	// https://stackoverflow.com/questions/9435385/split-a-string-using-c11
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}


set<string> load_barcodes(const string barcodes_path) {
	set<string> barcodes_set;
	string selected_barcodes_table = get_file_contents(barcodes_path.c_str());
	istringstream f(selected_barcodes_table);
	string line;
	int count = 0;
	while (getline(f, line)) {
		if (count>0) {
			vector<string> strs = split(line, ',');
			string firstelem = strs.at(0);
			barcodes_set.insert(firstelem);
		}
		count++;
	}
	return barcodes_set;
}


set<string> load_selected_read_names(const string read_names_path) {
	set<string> read_names_set;
	string read_names_string = get_file_contents(read_names_path.c_str());
	istringstream f(read_names_string);
	string line;
	while (getline(f, line)) {
		read_names_set.insert(line);
	}
	return read_names_set;
}



int main(int argc, char **argv) {
	if (argc == 4) {
		string read_names_path(argv[1]);
		string fasta_path(argv[2]);
		string barcodes_path(argv[3]);

		set<string> selected_read_names = load_selected_read_names(read_names_path);
		set<string> selected_barcodes = load_barcodes(barcodes_path);

		bool read_found = false;
		string read_name;


		FILE* fp = fopen(fasta_path.c_str(), "r");
		if (fp == NULL) {
			exit(EXIT_FAILURE);
		}

		char* line = NULL;
		size_t len = 0;
		int counter = 0;
		while ((getline(&line, &len, fp)) != -1) {
			if (counter % 2 == 0) {
				vector<string> split_line = split(line, ' ');
				read_name = split_line.at(0);
				read_name.erase(0, 1);
				if (selected_read_names.find(read_name)!=selected_read_names.end()) {
					read_found = true;
				} else {
					read_found = false;
				}

			} else if ((counter % 2 == 1) & (read_found == true)) {
				string barcode = string(line);
				barcode = barcode.substr(0, 16);
				if (selected_barcodes.find(barcode)!=selected_barcodes.end()) {
					cout << barcode << "," << read_name << endl;
				} else {
					cout << "NA" << "," << read_name << endl;
				}
			}
		    counter++;
		}
		fclose(fp);
		if (line) {
			free(line);
		}
		return 0;
	} else {
		cout << "C++ program for extracting Chromium barcode sequences of Chromium reads that have been selected for assembly" << endl;
		cout << "Wrong number of arguments" << endl;
		cout << "Usage:" << endl;
		cout << "Argument1: path to a text file with read names" << endl;
		cout << "Argument2: path to reads FASTA file, derived from FASTQ file (#1 file of paired end reads)" << endl;
		cout << "Argument3: path to a CSV file, with the first column containing sequences of barcodes that have more than one read associated with them. The CSV header line will be skipped by this program" << endl;
	}
}

