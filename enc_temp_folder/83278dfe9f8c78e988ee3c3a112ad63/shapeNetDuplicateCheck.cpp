// shapeNetDuplicateCheck.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "../PLY_Loader-master/plyloader.h"
#include "../libicp/src/icpPointToPlane.h"
#include "../libicp/src/icpPointToPoint.h"
#include "../point-to-plane/dist2Triangle.h"
#include "shapeNetDuplicateCheck.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <algorithm>
#include <mex.h>

using namespace std;

struct shapeNet
{
	string oriShapeId;
	string shapeId;
	double dist;
};

void createSamplePly(string datasetDir, string categoryDir);
void createRecomputeVNPly(string datasetDir, string categoryDir);
void encodeHash(string fileDir);
int findCommon(string s1, string s2);
void findKNN(string fileDir, int K);
void loadPlyMatrix(PLYModel plyModel, int plyNum, double* plyMatrix);
void icpMatching(string choice, mxArray* mx, int m, const Matrix &R, const Matrix &t);
void findDuplicate(string fileDir, string fileDir_VN, double thd);
int countDuplicate(string fileDir, double thd);
void loadPlyMxArray(PLYModel plyModel, string choice, int n, int m, mxArray* mx);
void loadPlyFaceNorm(mxArray* v, mxArray* f, int n, int m, mxArray* mx);
double maximum(double a, double b, double c);

void createSamplePly(string datasetDir, string categoryDir)
{
	ifstream infile;
	ofstream outfile;
	infile.open(datasetDir + categoryDir + "list.txt");

	string categoryName = "1000-" + categoryDir;
	categoryName.pop_back();
	string outfileName = datasetDir + "../createSamplePly_" + categoryName + ".bat";
	cout << outfileName << endl;
	outfile.open(outfileName);

	string fileName;

	int i = 0;
	while (getline(infile, fileName))
	{
		cout << "Processing #" << ++i << " file " << fileName << endl;
		outfile << "meshlabserver -i ShapeNetCore.v2/" << categoryDir << fileName << " -o ShapeNetCore.sample/1000-" << categoryDir << fileName << " -m vn -s poison_sample_mlx_VN.mlx" << endl;
	}

	infile.close();
	outfile.close();
}

void createRecomputeVNPly(string datasetDir, string categoryDir)
{
	ifstream infile;
	ofstream outfile;
	infile.open(datasetDir + categoryDir + "list.txt");

	string categoryName = categoryDir;
	categoryName.pop_back();
	string outfileName = datasetDir + "../createRecomputeVNPly_" + categoryName + ".bat";
	cout << outfileName << endl;
	outfile.open(outfileName);

	string fileName;

	int i = 0;
	while (getline(infile, fileName))
	{
		cout << "Processing #" << ++i << " file " << fileName << endl;
		outfile << "meshlabserver -i ShapeNetCore.v2/" << categoryDir << fileName << " -o ShapeNetCore.withVN/" << categoryDir << fileName << " -m vn -s recompute_norm_mlx.mlx" << endl;
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

void icpMatching(string choice, mxArray* mx, int m, const Matrix &R, const Matrix &t)
{
	double* ptr = mxGetPr(mx);
	double r00 = R.val[0][0]; double r01 = R.val[0][1]; double r02 = R.val[0][2];
	double r10 = R.val[1][0]; double r11 = R.val[1][1]; double r12 = R.val[1][2];
	double r20 = R.val[2][0]; double r21 = R.val[2][1]; double r22 = R.val[2][2];
	double t0 = t.val[0][0]; double t1 = t.val[1][0]; double t2 = t.val[2][0];

	if (choice.compare("v") == 0)
	{
		for (int i = 0; i < m; i++)
		{
			double sx = r00*ptr[i * 3 + 0] + r01*ptr[i * 3 + 1] + r02*ptr[i * 3 + 2] + t0; 
			double sy = r10*ptr[i * 3 + 0] + r11*ptr[i * 3 + 1] + r12*ptr[i * 3 + 2] + t1; 
			double sz = r20*ptr[i * 3 + 0] + r21*ptr[i * 3 + 1] + r22*ptr[i * 3 + 2] + t2; 
			ptr[i * 3 + 0] = sx;
			ptr[i * 3 + 1] = sy;
			ptr[i * 3 + 2] = sz;
		}
	}
	else if (choice.compare("n") == 0)
	{
		// G = (R^-1)^T for normal
		double det = r00*(r11*r22 - r12*r21) - r01*(r10*r22 - r12*r20) + r02*(r10*r21 - r11*r20);
		double m00 = (r11*r22 - r12*r21) / det;
		double m01 = -(r01*r22 - r02*r21) / det;
		double m02 = (r01*r12 - r02*r11) / det;
		double m10 = -(r10*r22 - r12*r20) / det;
		double m11 = (r00*r22 - r02*r20) / det;
		double m12 = -(r00*r12 - r02*r10) / det;
		double m20 = (r10*r21 - r11*r20) / det;
		double m21 = -(r00*r21 - r01*r20) / det;
		double m22 = (r00*r11 - r01*r10) / det;

		double m01t = m10; double m10t = m01; double m02t = m20; double m20t = m02; double m12t = m21; double m21t = m12;

		for (int i = 0; i < m; i++)
		{
			double sx = m00*ptr[i * 3 + 0] + m01t*ptr[i * 3 + 1] + m02t*ptr[i * 3 + 2];
			double sy = m10t*ptr[i * 3 + 0] + m11*ptr[i * 3 + 1] + m12t*ptr[i * 3 + 2];
			double sz = m20t*ptr[i * 3 + 0] + m21t*ptr[i * 3 + 1] + m22*ptr[i * 3 + 2];
			ptr[i * 3 + 0] = sx;
			ptr[i * 3 + 1] = sy;
			ptr[i * 3 + 2] = sz;
		}
	}
	else
	{
		cout << "Error usage of func icpMatching().\n";
		system("pause");
		exit(0);
	}
}

void findDuplicate(string fileDir, string fileDir_VN, double thd)
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
		
		string str = fileDir_VN + fileName;
		const char* plyName = str.c_str();
		cout << "Processing #" << iShape << " file " << fileName << endl;
		
		
		PLYModel plyModel(plyName, 1, 0);
		int plyNum = plyModel.vertexCount;
		int plyFaceNum = plyModel.faceCount;
		double plySize = maximum(plyModel.bvWidth, plyModel.bvHeight, plyModel.bvDepth);

		string str0 = fileDir + fileName;
		const char* plyName0 = str0.c_str();
		PLYModel plyModel0(plyName0, 1, 0);
		int plyNum0 = plyModel0.vertexCount;
		double* plyMatrix0 = (double*)calloc(3 * plyNum0, sizeof(double));
		loadPlyMatrix(plyModel0, plyNum0, plyMatrix0);
		IcpPointToPlane icp(plyMatrix0, plyNum0, 3);

		// point-to-plane parameter for 1st ply file
		mxArray *Ps = mxCreateDoubleMatrix(3, plyNum, mxREAL);
		loadPlyMxArray(plyModel, "v", 3, plyNum, Ps);
		mxArray *F = mxCreateDoubleMatrix(3, plyFaceNum, mxREAL);
		loadPlyMxArray(plyModel, "f", 3, plyFaceNum, F);
		mxArray *Nf = mxCreateDoubleMatrix(3, plyFaceNum, mxREAL);
		loadPlyFaceNorm(Ps, F, 3, plyFaceNum, Nf);
		mxArray *Nv = mxCreateDoubleMatrix(3, plyNum, mxREAL);
		loadPlyMxArray(plyModel, "n", 3, plyNum, Nv);
		mxArray *NormalCheck = mxCreateDoubleMatrix(1, 1, mxREAL);
		double* ptrNormalCheck = mxGetPr(NormalCheck);
		ptrNormalCheck[0] = 3.14;

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
				cout << " ="<< iShape <<"=> Neighbor #" << iCmpshape << " file " << cmpFilename << ": ";
				// compare with cmpFilename ply
				string str2 = fileDir + cmpFilename; 
				const char* plyName2 = str2.c_str();
				PLYModel plyModel2(plyName2, 1, 0);
				int plyNum2 = plyModel2.vertexCount;
				double* plyMatrix2 = (double*)calloc(3 * plyNum2, sizeof(double));
				loadPlyMatrix(plyModel2, plyNum2, plyMatrix2);

				Matrix R = Matrix::eye(3);
				Matrix t(3, 1);

				icp.fit(plyMatrix2, plyNum2, R, t, -1);

				// === icp package point-to-plane distance method =====================================================

				/*
				double* dist = new double(0);
				vector<int> inlierIdx = icp.getInliers(plyMatrix2, plyNum2, R, t, 10, dist);
				//cout << "       Distance: " << *dist << " , # of points: " << inlierIdx.back() << endl;
				shapeNet s;
				s.oriShapeId = fileName;
				s.shapeId = cmpFilename;
				s.dist = abs(*dist) / plyNum2;
				
				cout << abs(*dist) << ", ";
				*/
				// === Vangelis point-to-plane distance method ========================================================

				
				mxArray *Qs = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
				loadPlyMxArray(plyModel2, "v", 3, plyNum2, Qs);
				icpMatching("v", Qs, plyNum2, R, t);
				mxArray *Nq = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
				loadPlyMxArray(plyModel2, "n", 3, plyNum2, Nq);
				icpMatching("n", Nq, plyNum2, R, t);
				
				const mxArray *prhs[] = { Ps, F, Nf, Nv, Qs, Nq, NormalCheck };

				mxArray *Ds = mxCreateDoubleMatrix(1, plyNum2, mxREAL);
				loadPlyMxArray(plyModel, "0", 1, plyNum2, Ds);
				mxArray *plhs[] = { Ds };

				mexFunction(1, plhs, 7, prhs);

				double* ptr = mxGetPr(plhs[0]);
				double sptr = 0;
				for (int i = 0; i < plyNum2; i++)
				{
					sptr += ptr[i];
				}
				sptr = sptr / double(plyNum2);
				double absDis = sptr;
				double ratioDis = sptr / plySize;
				//cout << "Abs dis: " << absDis << " , ratio dis: " << ratioDis << endl;
				
				shapeNet s;
				s.oriShapeId = fileName;
				s.shapeId = cmpFilename;
				s.dist = ratioDis;

				
				// ====================================================================================================
				
				cout << s.dist << endl;

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

		free(plyMatrix0);
		cmpInfile.close();

		if (iShape % saveN == (saveN-1))
		{
			duplicatefile.close();
			duplicatefile_thd.close();
		}

		iShape++;
		//std::system("pause");
	}
	infile.close();
	knnfile.close();
	duplicatefile.close();
	duplicatefile_thd.close();
}


void findDuplicate_biDistance(string fileDir, string fileDir_VN, double thd)
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

		string str = fileDir_VN + fileName;
		const char* plyName = str.c_str();
		cout << "Processing #" << iShape << " file " << fileName << endl;


		PLYModel plyModel(plyName, 1, 0);   // mesh shape 1st ply -> for 1) point-to-plane dis
		int plyNum = plyModel.vertexCount;
		int plyFaceNum = plyModel.faceCount;
		double plySize = maximum(plyModel.bvWidth, plyModel.bvHeight, plyModel.bvDepth);

		string str0 = fileDir + fileName;
		const char* plyName0 = str0.c_str();
		PLYModel plyModel0(plyName0, 1, 0);   // sample point cloud 1st ply -> for 1) icp matching 2) inv-point-to-plane dis
		int plyNum0 = plyModel0.vertexCount;
		double* plyMatrix0 = (double*)calloc(3 * plyNum0, sizeof(double));
		loadPlyMatrix(plyModel0, plyNum0, plyMatrix0);
		IcpPointToPlane icp(plyMatrix0, plyNum0, 3);

		// point-to-plane parameter for 1st ply file
		mxArray *Ps = mxCreateDoubleMatrix(3, plyNum, mxREAL);
		loadPlyMxArray(plyModel, "v", 3, plyNum, Ps);
		mxArray *F = mxCreateDoubleMatrix(3, plyFaceNum, mxREAL);
		loadPlyMxArray(plyModel, "f", 3, plyFaceNum, F);


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

			if (iCmpshape != knnIdxint)
			{
				iCmpshape++;
				continue;
			}
			else
			{
				cout << " =" << iShape << "=> Neighbor #" << iCmpshape << " file " << cmpFilename << ": ";
				// compare with cmpFilename ply
				string str20 = fileDir + cmpFilename;
				const char* plyName20 = str20.c_str();
				PLYModel plyModel20(plyName20, 1, 0);   // sample point cloud 2nd ply -> for 1) icp matching 2) point-to-plane dis
				int plyNum20 = plyModel20.vertexCount;
				double* plyMatrix20 = (double*)calloc(3 * plyNum20, sizeof(double));
				loadPlyMatrix(plyModel20, plyNum20, plyMatrix20);

				Matrix R = Matrix::eye(3);
				Matrix t(3, 1);


				// ========================== point-to-plane dis ==============================
				// icp matching
				icp.fit(plyMatrix20, plyNum20, R, t, -1);  

				mxArray *Qs = mxCreateDoubleMatrix(3, plyNum20, mxREAL);
				loadPlyMxArray(plyModel20, "v", 3, plyNum20, Qs);
				icpMatching("v", Qs, plyNum20, R, t);

				const mxArray *prhs[] = { Ps, F, Qs };
				mxArray *Ds = mxCreateDoubleMatrix(1, plyNum20, mxREAL);
				loadPlyMxArray(plyModel20, "0", 1, plyNum20, Ds);
				mxArray *plhs[] = { Ds };

				// point-to-plane dis
				mexFunction_simple(1, plhs, 7, prhs);

				double* ptr = mxGetPr(plhs[0]);
				double sptr = 0;
				for (int i = 0; i < plyNum20; i++)
				{
					sptr += ptr[i];
				}
				sptr = sptr / double(plyNum20);
				double absDis = sptr;
				double ratioDis = sptr / plySize;
				//cout << "Abs dis: " << absDis << " , ratio dis: " << ratioDis << endl;

				shapeNet s;
				s.oriShapeId = fileName;
				s.shapeId = cmpFilename;
				s.dist = ratioDis;

				// ========================== inv-point-to-plane dis ==============================
				// icp matching
				IcpPointToPlane icp2(plyMatrix20, plyNum20, 3);
				icp2.fit(plyMatrix0, plyNum0, R, t, -1);

				// point-to-plane parameter for 2st ply file
				string str2 = fileDir_VN + cmpFilename;
				const char* plyName2 = str2.c_str();
				PLYModel plyModel2(plyName2, 1, 0);   // mesh shape 2nd ply -> for 1) inv-point-to-plane dis
				int plyNum2 = plyModel2.vertexCount;
				int plyFaceNum2 = plyModel2.faceCount;
				double plySize2 = maximum(plyModel2.bvWidth, plyModel2.bvHeight, plyModel2.bvDepth);
				
				mxArray *Ps_inv = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
				loadPlyMxArray(plyModel2, "v", 3, plyNum2, Ps_inv);
				mxArray *F_inv = mxCreateDoubleMatrix(3, plyFaceNum2, mxREAL);
				loadPlyMxArray(plyModel2, "f", 3, plyFaceNum2, F_inv);

				// query 1st sample points
				mxArray *Qs_inv = mxCreateDoubleMatrix(3, plyNum0, mxREAL);
				loadPlyMxArray(plyModel0, "v", 3, plyNum0, Qs_inv);
				icpMatching("v", Qs_inv, plyNum0, R, t);

				const mxArray *prhs_inv[] = { Ps_inv, F_inv, Qs_inv };
				mxArray *Ds_inv = mxCreateDoubleMatrix(1, plyNum0, mxREAL);
				loadPlyMxArray(plyModel0, "0", 1, plyNum0, Ds_inv);
				mxArray *plhs_inv[] = { Ds_inv };

				// point-to-plane dis
				mexFunction_simple(1, plhs_inv, 7, prhs_inv);

				ptr = mxGetPr(plhs_inv[0]);
				sptr = 0;
				for (int i = 0; i < plyNum0; i++)
				{
					sptr += ptr[i];
				}
				sptr = sptr / double(plyNum0);
				absDis = sptr;
				ratioDis = sptr / plySize2;
				//cout << "Abs dis: " << absDis << " , ratio dis: " << ratioDis << endl;

				s.dist = ratioDis;


				// =========================== END OF DISTANCE ==============================================================

				cout << s.dist << endl;

				cmpShapeNet.push_back(s);

				// free memory

				free(plyMatrix20);

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

		free(plyMatrix0);
		cmpInfile.close();

		if (iShape % saveN == (saveN - 1))
		{
			duplicatefile.close();
			duplicatefile_thd.close();
		}

		iShape++;
		//std::system("pause");
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
	else if (choice.compare("0") == 0)
	{
		for (int j = 0; j < m; j++)
		{
			ptr[j] = 0;
		}
	}
	else
	{
		cout << "Error usage of func loadPlyMxArray().\n";
		system("pause");
		exit(0);
	}
}

void loadPlyFaceNorm(mxArray* v, mxArray* f, int n, int m, mxArray* mx)
{
	// n = 3
	double* ptrv = mxGetPr(v);
	double* ptrf = mxGetPr(f);
	double* ptr = mxGetPr(mx);

	double v01[3], v02[3];
	int i0, i1, i2;

	for (int j = 0; j < m; j++)
	{
		i0 = ptrf[(3 * j) + 0];
		i1 = ptrf[(3 * j) + 1];
		i2 = ptrf[(3 * j) + 2];
		v01[0] = ptrv[(3 * i1) + 0] - ptrv[(3 * i0) + 0];
		v01[1] = ptrv[(3 * i1) + 1] - ptrv[(3 * i0) + 1];
		v01[2] = ptrv[(3 * i1) + 2] - ptrv[(3 * i0) + 2];
		v02[0] = ptrv[(3 * i2) + 0] - ptrv[(3 * i0) + 0];
		v02[1] = ptrv[(3 * i2) + 1] - ptrv[(3 * i0) + 1];
		v02[2] = ptrv[(3 * i2) + 2] - ptrv[(3 * i0) + 2];
		double* tmp = vec3cross(v01, v02);
		vec3normalize(tmp);
		for (int i = 0; i < 3; i++)
		{
			ptr[(3 * j) + i] = tmp[i];
		}
	}

}

double maximum(double a, double b, double c)
{
	double max = (a < b) ? b : a;
	return ((max < c) ? c : max);
}

int main(int argc, char** argv)
{
	std::cout << "Have " << argc << " arguments." << std::endl;
	if (argc < 3)
	{
		cout << "Usage: shapeNetDuplicateCheck.exe categoryNum choiceNum additionalParams.\n";
		cout << "Choice 1: Create poison-sample bash file & re-computer vertex norm bash file. \n";
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
	string param, param2;
	if (argc == 4)
		param = argv[3];
	if (argc == 5)
		param2 = argv[4];


	string datasetDir = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.v2/";
	string sampleDir = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.sample/";
	string datasetDir_VN = "E:/Dataset/Lun/Models/ShapeNetCore.v2/ShapeNetCore.withVN/";
	string categoryDir = categoryNum + "/";
	string fileDir = sampleDir + categoryDir;
	string fileDir_VN = datasetDir_VN + categoryDir;

	// Create poison sample bash file
	if (choiceNum.compare("1") == 0)
	{
		cout << "Creating poison sample (with vertex norm) bash file." << endl;
		createSamplePly(datasetDir, categoryDir);
		cout << "Creating RE-computer vertex norm bash file." << endl;
		createRecomputeVNPly(datasetDir, categoryDir);
	}
	if (choiceNum.compare("2") == 0)
	{
		cout << "Encode 3D point clouds to 6x6x6 grid representation." << endl;
		encodeHash(fileDir);
	}
	if (choiceNum.compare("3") == 0)
	{
		cout << "Casually find K most similar neighbors in the same category." << endl;
		findKNN(fileDir, 100);
	}
	if (choiceNum.compare("4") == 0)
	{
		cout << "Accuratly find distance between neighbors and judge duplicates." << endl;
		//findDuplicate(fileDir, fileDir_VN, atof(param.c_str()));
		findDuplicate_biDistance(fileDir, fileDir_VN, atof(param.c_str()));
	}
	
	if (choiceNum.compare("-2") == 0)
	{
		cout << "Testing tmp codes ..." << endl;
		mxArray *Ps = mxCreateDoubleMatrix(3, 3, mxREAL);
		double* ptrPs = mxGetPr(Ps);
		ptrPs[0] = 0; ptrPs[1] = 0; ptrPs[2] = 0;
		ptrPs[3 + 0] = 1; ptrPs[3 + 1] = 0; ptrPs[3 + 2] = 0;
		ptrPs[6 + 0] = 0; ptrPs[6 + 1] = 1; ptrPs[6 + 2] = 0;
		mxArray *F = mxCreateDoubleMatrix(3, 1, mxREAL);
		double* ptrF = mxGetPr(F);
		ptrF[0] = 0; ptrF[1] = 1; ptrF[2] = 2;
		mxArray *Nf = mxCreateDoubleMatrix(3, 1, mxREAL);
		double* ptrNf = mxGetPr(Nf);
		ptrNf[0] = 0; ptrNf[1] = 0; ptrNf[2] = 1;
		mxArray *Nv = mxCreateDoubleMatrix(3, 3, mxREAL);
		double* ptrNv = mxGetPr(Nv);
		ptrNv[0] = 0; ptrNv[1] = 0.3; ptrNv[2] = 0.9;
		ptrNv[3 + 0] = 0.4; ptrNv[3 + 1] = 0; ptrNv[3 + 2] = 0.9;
		ptrNv[6 + 0] = 0; ptrNv[6 + 1] = 0.2; ptrNv[6 + 2] = 0.9;
		mxArray *Qs = mxCreateDoubleMatrix(3, 1, mxREAL);
		double* ptrQs = mxGetPr(Qs);
		ptrQs[0] = -0.5; ptrQs[1] = 3; ptrQs[2] = -1;
		mxArray *Nq = mxCreateDoubleMatrix(3, 1, mxREAL);
		double* ptrNq = mxGetPr(Nq);
		ptrNq[0] = 0.4; ptrNq[1] = 0.3; ptrNq[2] = -0.9;
		mxArray *NormalCheck = mxCreateDoubleMatrix(1, 1, mxREAL);
		double* ptrNormalCheck = mxGetPr(NormalCheck);
		ptrNormalCheck[0] = 3.14;

		const mxArray *prhs[] = { Ps, F, Nf, Nv, Qs, Nq, NormalCheck };

		mxArray *Ds = mxCreateDoubleMatrix(1, 1, mxREAL);
		mxArray *plhs[] = { Ds };

		mexFunction(1, plhs, 7, prhs);

		double* ptr = mxGetPr(plhs[0]);
		cout << ptr[0] << endl;
	}


	if (choiceNum.compare("-1") == 0)
	{
		cout << "Testing tmp codes ..." << endl;

		string str0 = sampleDir + categoryDir + "103c9e43cdf6501c62b600da24e0965.ply";
		string str1 = datasetDir_VN + categoryDir + "103c9e43cdf6501c62b600da24e0965.ply";
		string str2 = sampleDir + categoryDir + "b4575e5e6161fd497b164268a44f7712.ply";
		const char* plyName0 = str0.c_str();
		const char* plyName1 = str1.c_str();
		const char* plyName2 = str2.c_str();
		PLYModel plyModel0(plyName0, 1, 0);
		PLYModel plyModel1(plyName1, 1, 0);
		PLYModel plyModel2(plyName2, 1, 0);
		int plyNum0 = plyModel0.vertexCount;
		int plyNum1 = plyModel1.vertexCount;
		int plyNum2 = plyModel2.vertexCount;
		int plyFaceNum1 = plyModel1.faceCount;

		double* plyMatrix0 = (double*)calloc(3 * plyNum0, sizeof(double));
		double* plyMatrix1 = (double*)calloc(3 * plyNum1, sizeof(double));
		double* plyMatrix2 = (double*)calloc(3 * plyNum2, sizeof(double));
		loadPlyMatrix(plyModel0, plyNum0, plyMatrix0);
		loadPlyMatrix(plyModel1, plyNum1, plyMatrix1);
		loadPlyMatrix(plyModel2, plyNum2, plyMatrix2);


		IcpPointToPlane icp(plyMatrix0, plyNum0, 3);
		Matrix R = Matrix::eye(3);
		Matrix t(3, 1);
		icp.fit(plyMatrix2, plyNum2, R, t, -1);

		cout << R.val[0][0] << " " << R.val[0][1] << " " << R.val[0][2] << endl;
		cout << R.val[1][0] << " " << R.val[1][1] << " " << R.val[1][2] << endl;
		cout << R.val[2][0] << " " << R.val[2][1] << " " << R.val[2][2] << endl;
		cout << t.val[0][0] << " " << t.val[0][1] << " " << t.val[0][2] << endl;

		mxArray *Ps = mxCreateDoubleMatrix(3, plyNum1, mxREAL);
		loadPlyMxArray(plyModel1, "v", 3, plyNum1, Ps);
		mxArray *F = mxCreateDoubleMatrix(3, plyFaceNum1, mxREAL);
		loadPlyMxArray(plyModel1, "f", 3, plyFaceNum1, F);
		mxArray *Nf = mxCreateDoubleMatrix(3, plyFaceNum1, mxREAL);
		loadPlyFaceNorm(Ps, F, 3, plyFaceNum1, Nf);
		mxArray *Nv = mxCreateDoubleMatrix(3, plyNum1, mxREAL);
		loadPlyMxArray(plyModel1, "n", 3, plyNum1, Nv);
		mxArray *Qs = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
		loadPlyMxArray(plyModel2, "v", 3, plyNum2, Qs);
		mxArray *Nq = mxCreateDoubleMatrix(3, plyNum2, mxREAL);
		loadPlyMxArray(plyModel2, "n", 3, plyNum2, Nq);
		mxArray *NormalCheck = mxCreateDoubleMatrix(1, 1, mxREAL);
		double* ptrNormalCheck = mxGetPr(NormalCheck);
		ptrNormalCheck[0] = 3.14;


		icpMatching("v", Qs, plyNum2, R, t);
		icpMatching("n", Nq, plyNum2, R, t);


		const mxArray *prhs[] = { Ps, F, Nf, Nv, Qs, Nq, NormalCheck };

		mxArray *Ds = mxCreateDoubleMatrix(1, plyNum2, mxREAL);
		loadPlyMxArray(plyModel1, "0", 1, plyNum2, Ds);
		mxArray *plhs[] = { Ds };

		mexFunction(1, plhs, 7, prhs);

		double* ptr = mxGetPr(plhs[0]);
		double sptr = 0;
		for (int i = 0; i < plyNum2; i++)
		{
			sptr += ptr[i];
		}
		sptr = sptr / double (plyNum2);
		double absDis = sptr;
		double ratioDis = sptr / maximum(plyModel1.bvWidth, plyModel1.bvHeight, plyModel1.bvDepth);
		cout << "Abs dis: " << absDis << " , ratio dis: " << ratioDis << endl;
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

