#include "ImageDataLoader.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

Eigen::MatrixXd ImageDataLoader::readImage(std::string path)
{
	ifstream myfile(path);

	Eigen::MatrixXd x;
	if (myfile.is_open())
	{
		char * memblock = new char[4];
		myfile.read(memblock, 4);
		
		std::stringstream str(memblock);
		str >> numRows;
		if (!str)
		{
			cout << "conv failed." << endl;
		}

		//cout << "rows = " << numRows << endl;

		myfile.read(memblock, 1); //delimmiter 

		myfile.read(memblock, 4);

		std::stringstream str2(memblock);
		str2 >> numCols;
		if (!str2)
		{
			cout << "conv failed." << endl;
		}

		//cout << "cols = " << numCols << endl; 

		myfile.read(memblock, 1); //delimmiter 


		Eigen::MatrixXd data(numRows, numCols);
		for (int i = 0; i < numRows; i++)
		{
			for (int j = 0; j < numCols; j++)
			{
				myfile.read(memblock, 1);

				std::stringstream str3(memblock);
				str3 >> data(i,j);
				if (!str2)
				{
					cout << "conv failed." << endl;
				}

			}
			cout << endl;
		}
		
		x = data;
		myfile.close();
	}
	else cout << "Unable to open file";

	//cout << x;

	return x;
}

std::vector<Eigen::MatrixXd> ImageDataLoader::readAllImages(std::string path, size_t numImages)
{
	std::vector<Eigen::MatrixXd> images(numImages);
	int imNum = 0;
	for (int i = 0; i < numImages; i++)
	{
		stringstream fileNameFormat;
		fileNameFormat << path <<"/slice" << setfill('0') << setw(4) << imNum+1 << ".txt";
		//fileNameFormat << path <<"/square" << setfill('0') << setw(2) << i+1 << ".txt";
		//fileNameFormat << path << "/lattice_test_" << setfill('0') << setw(2) << i+1 << ".txt";
		string fileName;
		fileNameFormat >> fileName;

		imagePath.push_back(fileName);
		images[i] = readImage(fileName);
		imNum+=3;
	}
	//manually add the last slice

			stringstream fileNameFormat;
	fileNameFormat << path << "/slice" << setfill('0') << setw(4) << 111 << ".txt";
	string fileName;
	fileNameFormat >> fileName;
	imagePath.push_back(fileName);
	images.push_back(readImage(fileName));

		//stringstream fileNameFormat;
	fileNameFormat << path << "/slice" << setfill('0') << setw(4) << 112 << ".txt";
	//string fileName;
	fileNameFormat >> fileName;
	imagePath.push_back(fileName);
	images.push_back(readImage(fileName));


	//stringstream fileNameFormat;
	fileNameFormat << path << "/slice" << setfill('0') << setw(4) << 113 << ".txt";
	//string fileName;
	fileNameFormat >> fileName;
	imagePath.push_back(fileName);
	images.push_back(readImage(fileName));
	return images;
}