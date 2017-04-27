// shapeNetDuplicateCheck.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "../PLY_Loader-master/plyloader.h"
#include "../libicp/src/icpPointToPlane.h"
#include "../libicp/src/icpPointToPoint.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include<algorithm>
#include <mex.h>

using namespace std;

struct shapeNet
{
	string oriShapeId;
	string shapeId;
	double dist;
	int inlierNum;
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
		outfile << "meshlabserver -i ShapeNetCore.v2/" << categoryDir << fileName << " -o ShapeNetCore.sample/" << categoryDir << fileName << " -m vn -s poison_sample_mlx.mlx" << endl;
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

void loadPlyMatrix(PLYModel plyModel, int plyNum, double* plyMatrix)
{
	for (int i = 0; i < plyNum; i++)
	{
		plyMatrix[i * 3 + 0] = plyModel.positions[i].x;
		plyMatrix[i * 3 + 1] = plyModel.positions[i].y;
		plyMatrix[i * 3 + 2] = plyModel.positions[i].z;
	}
}

void findDuplicate(string fileDir, double thd)
{
	int startIdx = 0;
	ifstream infile, cmpInfile, knnfile;
	infile.open(fileDir + "list.txt");
	knnfile.open(fileDir + "knn.txt");

	ofstream duplicatefile, duplicatefile_thd;
	int saveN = 100;
	string fileName;
	string knnString;

	int iShape = 0;
	while (getline(infile, fileName))
	{
		// read knn file
		if (!getline(knnfile, knnString))
		{
			cout << "EOF in knn.txt file!\n";
			break;
		}

		// start iShape index
		if (iShape < startIdx)
		{
			cout << "Skip # " << iShape << "file " << fileName << endl;
			iShape++;
			continue;
		}

		// read 1st ply file
		string str = fileDir + fileName;
		const char* plyName = str.c_str();
		cout << "Processing #" << iShape << " file " << fileName << endl;
		

		PLYModel plyModel(plyName, 0, 0);
		int plyNum = plyModel.vertexCount;
		double* plyMatrix = (double*)calloc(3 * plyNum, sizeof(double));
		loadPlyMatrix(plyModel, plyNum, plyMatrix);
		IcpPointToPlane icp(plyMatrix, plyNum, 3);

		// process knn file
		istringstream knnStringstream(knnString);
		string knnIdx;
		knnStringstream >> knnIdx; // exclude knn distance
		knnStringstream >> knnIdx; // read first neighbor
		int knnIdxint;
		knnIdxint = stoi(knnIdx);

		// ofstream duplicate file
		if (iShape % saveN == 0)
		{
			duplicatefile.open(fileDir + "dis/duplicate-" + to_string(iShape) + ".txt");
			duplicatefile_thd.open(fileDir + "dis/duplicate_thd-" + to_string(iShape) + ".txt");
		}
		duplicatefile << "ori " << iShape << " " << fileName << endl;
		duplicatefile_thd << "ori " << iShape << " " << fileName << endl;

		// build distance - neighborId sorting pair
		vector<shapeNet> cmpShapeNet;

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
				cout << " ="<< iShape <<"=> Neighbor #" << iCmpshape << " file " << cmpFilename << ": \n";
				// compare with cmpFilename ply
				string str2 = fileDir + cmpFilename; 
				const char* plyName2 = str2.c_str();
				PLYModel plyModel2(plyName2, 0, 0);
				int plyNum2 = plyModel2.vertexCount;
				//int plyNum2 = min(plyNum, plyModel2.vertexCount);
				double* plyMatrix2 = (double*)calloc(3 * plyNum2, sizeof(double));
				loadPlyMatrix(plyModel2, plyNum2, plyMatrix2);

				Matrix R = Matrix::eye(3);
				Matrix t(3, 1);

				icp.fit(plyMatrix2, plyNum2, R, t, -1);

				// === icp package point-to-plane distance method ===

				//double* dist = new double(0);
				//vector<int> inlierIdx = icp.getInliers(plyMatrix2, plyNum2, R, t, 10, dist);
				//cout << "       Distance: " << *dist << " , # of points: " << inlierIdx.back() << endl;

				// === Vangelis point-to-plane distance method ===
				mxArray **plhs;
				mxArray **prhs;
				//mexFunction(1, plhs, 7, prhs);

				shapeNet s;
				s.oriShapeId = fileName;
				s.shapeId = cmpFilename;
				//s.dist = abs(*dist);
				//s.inlierNum = inlierIdx.back();
				cmpShapeNet.push_back(s);

				// free memory
				
				free(plyMatrix2);

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


				cout << "Test for only 1\n";
				break;

			}
		}

		// sort cmpShapeNet according to dist
		std::sort(cmpShapeNet.begin(), cmpShapeNet.end(), [](const auto& i, const auto& j) { return i.dist < j.dist; });

		cout << " ====> For original file " << fileName << " after sorting: \n";
		for (int i = 0; i < cmpShapeNet.size(); i++)
		{
			cout << "       Distance: " << cmpShapeNet[i].dist << " , Shape Id: " << cmpShapeNet[i].shapeId << endl;
			duplicatefile << "cmp " << cmpShapeNet[i].dist << " " << cmpShapeNet[i].shapeId << endl;
			if (cmpShapeNet[i].dist < thd)
			{
				duplicatefile_thd << "cmp " << cmpShapeNet[i].dist << " " << cmpShapeNet[i].shapeId << endl;
			}
		}

		free(plyMatrix);
		cmpInfile.close();

		if (iShape % saveN == (saveN-1))
		{
			duplicatefile.close();
			duplicatefile_thd.close();
		}

		iShape++;
		//system("pause");
	}
	infile.close();
	knnfile.close();
	duplicatefile.close();
	duplicatefile_thd.close();
}

int countDuplicate(string fileDir, double thd)
{
	ifstream infile, duplicatefile;
	infile.open(fileDir + "dis/list.txt");
	ofstream outfile;
	outfile.open(fileDir + "dis/count.txt");
	string fileName;

	int totalDuplicate = 0;
	int fileDuplicate = 0;
	bool f = false;
	while (getline(infile, fileName))
	{
		string lines;
		duplicatefile.open(fileDir + "dis/" + fileName);
		while (getline(duplicatefile, lines))
		{
			if (lines[0] == 'o')
			{
				if (f)
				{
					outfile << fileDuplicate << endl;
					cout << fileDuplicate << endl;
				}
				istringstream lineStream(lines);
				string parts;
				lineStream >> parts;
				lineStream >> parts;
				outfile << parts << ' ';
				cout << parts << ' ';
				lineStream >> parts;
				outfile << parts << ' ';
				cout << parts << ' ';
				fileDuplicate = 0;
				
			}
			if (lines[0] == 'c')
			{
				istringstream lineStream(lines);
				string parts;
				lineStream >> parts;
				lineStream >> parts;
				double dis = atof(parts.c_str());
				if (dis < thd)
				{
					fileDuplicate++;
					totalDuplicate++;
				}
				f = true;
			}
			
		}


		duplicatefile.close();
	}
	outfile << fileDuplicate << endl;
	cout << fileDuplicate << endl;
	outfile << "Total duplicate pairs: " << totalDuplicate << endl;
	infile.close();
	outfile.close();
	return totalDuplicate;
}

void loadPlyMxArray(PLYModel plyModel, string choice, int n, int m, mxArray* mx)
{
	// n = 3

	double* ptr = mxGetPr(mx);
	if (choice.compare("v") == 0)
	{
		for (int j = 0; j < m; j++)
		{
			ptr[(3 * j) + 0] = plyModel.positions[j].x;
			ptr[(3 * j) + 1] = plyModel.positions[j].y;
			ptr[(3 * j) + 2] = plyModel.positions[j].z;
		}
	}
	else if (choice.compare("n") == 0)
	{
		for (int j = 0; j < m; j++)
		{
			ptr[(3 * j) + 0] = plyModel.normals[j].x;
			ptr[(3 * j) + 1] = plyModel.normals[j].y;
			ptr[(3 * j) + 2] = plyModel.normals[j].z;
		}
	}
	else if (choice.compare("f") == 0)
	{
		for (int j = 0; j < m; j++)
		{
			ptr[(3 * j) + 0] = plyModel.faces[j].x;
			ptr[(3 * j) + 1] = plyModel.faces[j].y;
			ptr[(3 * j) + 2] = plyModel.faces[j].z;
		}
	}
	else
	{
		cout << "Error usage of func loadPlyMxArray().\n";
		system("pause");
		exit(0);
	}
}


int main(int argc, char** argv)
{
	std::cout << "Have " << argc << " arguments." << std::endl;
	if (argc < 3)
	{
		cout << "Usage: shapeNetDuplicateCheck.exe categoryNum choiceNum additionalParams.\n";
		cout << "Choice 1: Create poison sample bash file. \n";
		cout << "Choice 2: Encode to 6x6x6 grid representation. \n";
		cout << "Choice 3: Casually find K most similar neighbors in the same category. \n";
		cout << "Choice 4: Accuratly find distance between neighbors and judge duplicates. \n";
		cout << "Choice 5: Total num of duplicates based on thd. \n";

		return 0;
	}
	for (int i = 0; i < argc; ++i) {
		std::cout << argv[i] << std::endl;
	}
	string categoryNum = argv[1];
	string choiceNum = argv[2];
	if (argc == 4)
		string param = argv[3];
	if (argc == 5)
		string param2 = argv[4];


	string datasetDir = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.v2/";
	string sampleDir = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.sample/";
	string categoryDir = categoryNum + "/";
	string fileDir = sampleDir + categoryDir;

	// Create poison sample bash file
	if (choiceNum.compare("1") == 0)
	{
		cout << "Creating poison sample bash file." << endl;
		createSamplePly(datasetDir, categoryDir);
	}
	if (choiceNum.compare("2") == 0)
	{
		cout << 22222 << endl;
	}
	if (choiceNum.compare("3") == 0)
	{
		cout << 3333 << endl;
	}
	if (choiceNum.compare("4") == 0)
	{
		cout << 44444 << endl;
	}
	

	if (choiceNum.compare("-1") == 0)
	{
		cout << "Testing tmp codes ..." << endl;

		string str1 = datasetDir + categoryDir + "10155655850468db78d106ce0a280f87.ply";
		string str2 = fileDir + "e6f0811f15286120cedbd07f4cf21a81.ply";
		const char* plyName1 = str1.c_str();
		const char* plyName2 = str2.c_str();
		PLYModel plyModel1(plyName1, 1, 0);
		PLYModel plyModel2(plyName2, 1, 0);
		int plyNum1 = plyModel1.vertexCount;
		int plyNum2 = plyModel2.vertexCount;
		int plyFaceNum1 = plyModel1.faceCount;
		cout << "p1: v=" << plyNum1 << " f=" << plyFaceNum1 << " , p2: v=" << plyNum2 << endl;
		double* plyMatrix1 = (double*)calloc(3 * plyNum1, sizeof(double));
		double* plyMatrix2 = (double*)calloc(3 * plyNum2, sizeof(double));
		loadPlyMatrix(plyModel1, plyNum1, plyMatrix1);
		loadPlyMatrix(plyModel2, plyNum2, plyMatrix2);

		for(int i=0; i<5; i++)
			cout << plyModel1.faces[0].x << "," << plyModel1.faces[0].y << "," << plyModel1.faces[0].z << endl;


		mxArray *Ps = mxCreateDoubleMatrix(3, plyNum1, mxREAL);
		loadPlyMxArray(plyModel1, "v", 3, plyNum1, Ps);
		mxArray *F = mxCreateDoubleMatrix(3, plyFaceNum1, mxREAL);
		loadPlyMxArray(plyModel1, "f", 3, plyFaceNum1, Ps);
		mxArray *Nf = mxCreateDoubleMatrix(3, plyFaceNum1, mxREAL);
		// =========================== Face norm !!!!!!!!======================
		mxArray *Nv = mxCreateDoubleMatrix(3, plyNum1, mxREAL);
		loadPlyMxArray(plyModel1, "n", 3, plyNum1, Ps);
		mxArray *Qs = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
		loadPlyMxArray(plyModel2, "v", 3, plyNum2, Ps);
		mxArray *Nq = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
		// =========================== vertex norm for sample points!!!!!!!!======================
		mxArray *NormalCheck = mxCreateDoubleMatrix(1, 1, mxREAL);
		// =========================== what is norm check??? !!!!!!!!======================
		const mxArray *prhs[] = { Ps, F, Nf, Nv, Qs, Nq, NormalCheck };


		mxArray *b0 = mxCreateDoubleMatrix(1, 100, mxREAL);
		mxArray *plhs[] = {b0};



		cout << "444" << endl;

		mexFunction(1, plhs, 7, prhs);
	}
	
	

	// Encode to 6x6x6 grid representation
	//encodeHash(fileDir);
	
	// Casually find K most similar neighbors in the same category.
	//findKNN(fileDir, 100);

	// Accuratly find distance between neighbors and judge duplicates.
	//findDuplicate(fileDir, 20);

	// Total num of duplicates based on thd.
	//int c = countDuplicate(fileDir, 20);
	//cout << c << endl;
	


	system("pause");


    return 0;
}

