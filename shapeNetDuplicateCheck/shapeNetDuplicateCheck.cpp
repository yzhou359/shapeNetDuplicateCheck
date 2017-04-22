// shapeNetDuplicateCheck.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "../PLY_Loader-master/plyloader.h"
#include "../libicp/src/icpPointToPlane.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include<algorithm>

using namespace std;

class shapeNet
{
public:
	string shapeId;
	string hashCode;
	int neighborDist;
	int* neighborIdxs;
};

void createSamplePly(string datasetDir, string categoryDir)
{
	ifstream infile;
	ofstream outfile;
	infile.open(datasetDir + categoryDir + "list.txt");

	string categoryName = categoryDir;
	categoryName.pop_back();
	string outfileName = datasetDir + "../createSamplePly_" + categoryName + ".bat";
	cout << outfileName << endl;
	outfile.open(outfileName);

	string fileName;

	int i = 0;
	while (getline(infile, fileName))
	{
		cout << "Processing #" << ++i << " file " << fileName << endl;
		outfile << "meshlabserver -i ShapeNetCore.v2/" << categoryDir << fileName << " -o ShapeNetCore.sample/" << categoryDir << fileName << " -s poison_sample_mlx.mlx" << endl;
	}

	infile.close();
	outfile.close();
}

void encodeHash(string fileDir)
{
	const int r = 8;
	char hashCode[r*r*r+1];

	ifstream infile;
	infile.open(fileDir + "list.txt");
	ofstream hashfile;
	hashfile.open(fileDir + "hash.txt");

	string fileName;

	int iShape = 0;
	while (getline(infile, fileName))
	{
		for (int i = 0; i < (r*r*r); i++) hashCode[i] = '0';
		hashCode[r*r*r] = '\0';

		string str = fileDir + fileName;
		const char* plyName = str.c_str();
		cout << "Processing #" << ++iShape << " file " << fileName << endl;
		PLYModel tmpPly(plyName, 0, 0);
		
		for (int i = 0; i < tmpPly.vertexCount; i++)
		{
			glm::vec3 tmpVec = tmpPly.positions[i] - tmpPly.min;
			glm::vec3 scaleVec = (tmpPly.max - tmpPly.min) + glm::vec3(0.000001,0.000001,0.000001);
			int idx = int(tmpVec.x / scaleVec.x*r)*r*r + int(tmpVec.y / scaleVec.y*r)*r + int(tmpVec.z / scaleVec.z*r);
			if (idx > r*r*r)
			{
				cout << idx << ": " << int(tmpVec.x / scaleVec.x*r) << "," << int(tmpVec.y / scaleVec.y*r) << "," << int(tmpVec.z / scaleVec.z*r) << endl;
				system("pause");
			}
			hashCode[idx] = '1';
			
		}
		hashfile << hashCode << endl;

		tmpPly.FreeMemory();
		

		//break;
	}

	infile.close();
	hashfile.close();
}

int findCommon(string s1, string s2)
{
	int count = 0;
	for (int i = 0; i < s1.length(); i++) 
	{
		if (s1[i] == s2[i]) 
		{
			count++;
		}
	}
	return count;
}

void findKNN(string fileDir, int K)
{
	ofstream knnfile;
	knnfile.open(fileDir + "knn.txt");

	ifstream hashfile, cmpHashfile;
	hashfile.open(fileDir + "hash.txt");
	string hashCode;
	int n_line = 0;
	while (getline(hashfile, hashCode))
	{
		n_line++;
	}
	hashfile.close();

	hashfile.open(fileDir + "hash.txt");
	int iShape = 0;
	while (getline(hashfile, hashCode))
	{
		cout << "Processing #" << ++iShape << endl;

		cmpHashfile.open(fileDir + "hash - Copy.txt");
		if (!cmpHashfile.is_open())
		{
			cout << "Cannot open cmpHashfile!" << endl;
			break;
		}

		string cmpHashcode;
		int* dist = new int[n_line];
		int* distCopy = new int[n_line];
		int cmpIdx = 0;
		while (getline(cmpHashfile, cmpHashcode))
		{
			int fc = findCommon(hashCode, cmpHashcode);
			dist[cmpIdx] = fc;
			distCopy[cmpIdx] = fc;
			cmpIdx++;
		}
		sort(dist, dist + n_line);

		int boundDist = dist[n_line-K];

		knnfile << boundDist << " ";

		for (int i = 0; i < n_line; i++)
		{
			if (distCopy[i] > boundDist)
			{
				knnfile << i << " ";
			}
		}
		knnfile << endl;

		cmpHashfile.close();
	}
	hashfile.close();
}

void findDuplicate(string fileDir, double thd)
{
	ifstream infile, cmpInfile, knnfile;
	infile.open(fileDir + "list.txt");
	knnfile.open(fileDir + "knn.txt");

	string fileName;
	string knnString;

	int iShape = 0;
	while (getline(infile, fileName))
	{
		// read 1st ply file
		string str = fileDir + fileName;
		const char* plyName = str.c_str();
		cout << "Processing #" << iShape << " file " << fileName << endl;
		PLYModel tmpPly(plyName, 0, 0);




		if (!getline(knnfile, knnString))
		{
			cout << "EOF in knn.txt file!\n";
			break;
		}
		istringstream knnStringstream(knnString);
		string knnIdx;
		knnStringstream >> knnIdx; // exclude knn distance
		knnStringstream >> knnIdx; // read first neighbor
		int knnIdxint;
		knnIdxint = stoi(knnIdx);

		// read 2nd ply file and judge whether it's neighbor
		int iCmpshape = 0;
		string cmpFilename;
		cmpInfile.open(fileDir + "list - Copy.txt");
		while (getline(cmpInfile, cmpFilename))
		{
			if (knnIdxint == iShape)
			{
				knnStringstream >> knnIdx; // read next neighbor
				try
				{
					knnIdxint = stoi(knnIdx);
				}
				catch (exception e)
				{
					cout << "End of knn Id.\n";
					break;
				}
			}
			
			if(iCmpshape != knnIdxint)
			{
				iCmpshape++;
				continue;
			}
			else
			{
				cout << iCmpshape << ": " << cmpFilename << endl;
				// compare with cmpFilename ply













				iCmpshape++;
				knnStringstream >> knnIdx; // read next neighbor
				try
				{
					knnIdxint = stoi(knnIdx);
				}
				catch (exception e)
				{
					cout << "End of knn Id.\n";
					break;
				}
			}
		}

		cmpInfile.close();
		iShape++;
		system("pause");


	}

}



int main()
{

	string datasetDir = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.v2/";
	string sampleDir = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.sample/";
	string categoryDir = "02691156/";
	string fileDir = sampleDir + categoryDir;
	
	// Create poison sample bash file
	//createSamplePly(datasetDir, categoryDir);

	// Encode to 6x6x6 grid representation
	//encodeHash(fileDir);
	
	// Casually find K most similar neighbors in the same category.
	//findKNN(fileDir, 100);

	// Accuratly find distance between neighbors and judge duplicates.
	//findDuplicate(fileDir, 0.1);

	int32_t dim = 3;
	int32_t num = 10000;

	// allocate model and template memory
	double* M = (double*)calloc(3 * num, sizeof(double));
	double* T = (double*)calloc(3 * num, sizeof(double));

	// set model and template points
	cout << endl << "Creating model with 10000 points ..." << endl;
	cout << "Creating template by shifting model by (1,1,1) ..." << endl;
	int32_t k = 0;
	for (double x = -2; x<2; x += 0.04) {
		for (double y = -2; y<2; y += 0.04) {
			double z = 5 * x*exp(-x*x - y*y);
			M[k * 3 + 0] = x;
			M[k * 3 + 1] = y;
			M[k * 3 + 2] = z;
			T[k * 3 + 0] = x - 1;
			T[k * 3 + 1] = y - 1;
			T[k * 3 + 2] = z - 1;
			k++;
		}
	}

	// start with identity as initial transformation
	// in practice you might want to use some kind of prediction here
	Matrix R = Matrix::eye(3);
	Matrix t(3, 1);

	// run point-to-plane ICP (-1 = no outlier threshold)
	cout << endl << "Running ICP (point-to-plane, no outliers)" << endl;
	IcpPointToPlane icp(M, num, dim);
	icp.fit(T, num, R, t, -1);

	// results
	cout << endl << "Transformation results:" << endl;
	cout << "R:" << endl << R << endl << endl;
	cout << "t:" << endl << t << endl << endl;

	// free memory
	free(M);
	free(T);


	system("pause");


    return 0;
}

