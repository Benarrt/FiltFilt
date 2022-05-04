#include <iostream>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#include "FiltFilt.h"


int main(int argc, char **argv) 
{
	if(argc == 1)
	{
		std::cout << "Usage:\n"
		<< "-f [Data file name]\n"
		<< "-d [Row delimeted]\n"
		<< "-c [Column index (from 0)]\n"
		<< "-r [Row to start from (indexed from 0)]\n"
		<< "e.g. -f data.txt -d , -c 1 -r 1\n"
		<< "See attached example data files.\nData with single colum (on record per line) can ignore -d and -c params"
		<< std::endl;
		return 0;
	}

	int column = 0;
	int row = 0;
	std::string delim;
	std::string filename; 
	std::vector<std::vector<double>> data;
	for(int i = 1; i < argc-1; i++)
	{
		std::string arg(argv[i]);
		if(arg == "-f")
		{
			filename = argv[++i];
		}
		else if(arg == "-d")
		{
			delim = argv[++i];
		}
		else if(arg == "-c")
		{
			column = std::stoi(argv[++i]);
		}
		else if(arg == "-r")
		{
			row = std::stoi(argv[++i]);
			if(row < 0)
				row = 0;
		}
	}

	std::string dataRow;
	std::ifstream inputFile;
	inputFile.open(filename);
	if(!inputFile.is_open())
	{
		std::cout << "File " << filename << " not found" << std::endl;
		return 0;
	}

	while (getline(inputFile, dataRow))
	{
		if(row > 0)
		{
			row--;
			continue;
		}
		int curColumn = 0;
		if(data.size() <= curColumn)
		{
			data.push_back({});
		}

		if(delim.empty())
		{
			data[curColumn].push_back(std::stod(dataRow));
		}
		else
		{
			int lastPost = 0;
			int pos = dataRow.find(delim);
			while(pos != std::string::npos)
			{
				if(data.size() <= curColumn)
				{
					data.push_back({});
				}

				data[curColumn].push_back(std::stod(dataRow.substr(lastPost, pos-lastPost)));
				lastPost = pos+1;
				pos = dataRow.find(delim, pos+1);
				curColumn++;
			}
		}
	}
	inputFile.close();

	if(column > data.size())
	{
		std::cout << "File " << filename << " has less than " << column+1 << " columns" << std::endl;
		return 0;
	}

	std::vector<double> b;
	std::vector<double> a{1};

	int res = system("python3 ./genCoeff.py");
	std::ifstream bCoeffFile;
	bCoeffFile.open("./bCoeff");
	if(!bCoeffFile.is_open())
	{
		std::cout << "File " << "./bCoeff" << " not found" << std::endl;
		return 0;
	}

	while (getline(bCoeffFile, dataRow))
	{
		b.push_back(std::stod(dataRow));
	}

	bCoeffFile.close();

	std::ifstream aCoeffFile;
	aCoeffFile.open("./aCoeff");
	if(!aCoeffFile.is_open())
	{
		std::cout << "File " << "./aCoeff" << " not found" << std::endl;
		return 0;
	}

	while (getline(aCoeffFile, dataRow))
	{
		a.push_back(std::stod(dataRow));
		std::cout << a.back() << std::endl;
	}

	aCoeffFile.close();

	std::vector<double> outY;
	Digital::FiltFilt sdfiltfilt;
	auto t_start = std::chrono::high_resolution_clock::now();
	sdfiltfilt(b, a, data[column], outY);
	auto t_end = std::chrono::high_resolution_clock::now();
	auto elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();

	std::cout << "C++ filtfilt took " << elapsed_time_ms << " ms" << std::endl;
	
	std::ofstream outYFile;
	outYFile.open("./demoFiltFilt");
	if(!outYFile.is_open())
	{
		std::cout << "File " << "./demoFiltFilt" << " not found" << std::endl;
		return 0;
	}

	bool first = true;
	for(auto y : outY)
	{
		if(!first)
		{
			outYFile << "\n";
		}
		else
		{
			first = !first;
		}
		
		outYFile << std::to_string(y);
	}

	outYFile.close();

	res = system(std::string("python3 ./filtfilt.py " + filename + " " + delim + " " + (column == 0 ? "" : std::to_string(column))).c_str());

	return res;
}
