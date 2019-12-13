/*
Program to remove barcodes and linkers from FASTQ files of Chromium reads.
Author: Eerik Aunin. Copyright (C) 2018, 2019 Genome Research Ltd.
This program is distributed under the terms of the GNU General Public License.
Argument1: path to FASTQ file of reads #1 of Chromium paired end reads.
Output: FASTQ file where the first 22 nucleotides (bar code and linker) have been deleted.

Compiling (requires g++ version >= 4.3, has been tested with gcc-8.2.0 on Ubuntu Linux):
g++ -std=c++11 chromium_reads_remove_barcodes.cpp -o chromium_reads_remove_barcodes
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
using namespace std;

int main(int argc,  char **argv) {
	if (argc == 2) {
		std::string file_path(argv[1]);
		std::ifstream file(file_path.c_str());
		if (file.good()) {
			std::string str;
			int counter = 0;
			int bar_code_and_linker_length = 22;

			while (std::getline(file, str)) {
				int size = str.length();
				int line_type = counter % 4;
				if (line_type == 1 || line_type == 3) {
					string str2;
					if (size > bar_code_and_linker_length) {
						str2 = str.substr(bar_code_and_linker_length + 1, size - bar_code_and_linker_length);
					} else {
						cout << "Encountered a sequence that is too short" << endl;
						break;
						exit(EXIT_FAILURE);
					}
					cout << str2 << endl;

				} else {
					cout << str << endl;
				}
				++counter;
			}
		} else {
			cout << file_path << ": file not found" << endl;
			exit(EXIT_FAILURE);
		}
	} else {
		cout << argv[0] << ": tool for removing barcodes and linkers from Chromium read FASTQ files" << endl;
		cout << "Usage: ";
		cout << argv[0] << " [FASTQ file]" << endl;
		exit(EXIT_FAILURE);
	}

    return 0;
}
