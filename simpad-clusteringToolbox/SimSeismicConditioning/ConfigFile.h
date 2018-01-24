#ifndef CONFIG_FILE_READER
#define CONFIG_FILE_READER

#include <string>
#include <time.h>
#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>

#include "UniversalSGEMSreader.h"
#include "FileReadingToolBox.h"

using namespace std;

class ConfigFile
{
public:
	ConfigFile(char* path);
	~ConfigFile();
	
	int GetTIWidth();
	int GetTIHeight();
	int GetTIDepth();

	int GetTemplateWidth();
	int GetTemplateHeight();
	int GetTemplateDepth();

	int GetRealizationWidth();
	int GetRealizationHeight();
	int GetRealizationDepth();

	int GetOverlapWidth();
	int GetOverlapHeight();
	int GetOverlapDepth();

	int GetAmountOfSimulations();
	double GetMinimumCompatibility();
	double GetSigmaSquared();

	short*** GetTrainingImage();
	double*** GetSyntheticSeismicData();
	double*** GetSeismicConditioning();

	double GetThresholdMaxDistance();
	int GetMaxAmountOfPatternsInTheQueue();
	int GetMaxCandidates();

	string GetRealizationFileName11();


	template <class T>
	static void GetSubMatrix(T*** bigMatrix, int x, int y, int z, T*** smallMatrix, int matrixWidth, int matrixHeight, int matrixDepth);
	static short*** InitializeMatrix_short(int width, int height, int depth);
	static double*** InitializeMatrix_double(int width, int height, int depth);
	static void DeleteMatrix(short*** simulationGrid, int width, int height, int depth);
	static void DeleteMatrix(double*** simulationGrid, int width, int height, int depth);


private:
	//char* path1;
	//string GetNextValidLine(ifstream& fs);
	
	
	int tiWidth  = 50;
	int tiHeight = 50;
	int tiDepth  = 50;

	int templateWidth  = 16;
	int templateHeight = 16;
	int templateDepth  = 16;

	int realizationWidth  = 76;//52;
	int realizationHeight = 76;//52;
	int realizationDepth  = 76;//52;

	int overlapWidth  = 4;
	int overlapHeight = 4;
	int overlapDepth  = 4;

	int amountOfSimulations = 3;
	double minimumCompatibility = 0.9;
	double sigmaSquared = 100;

	//string tiPath = "TI_3D\\fold_categorical_SMALL.SGEMS";
	//string seismicPath = "TI_3D\\fold_convoluted32_SMALL.sgems";
	//string syntheticPath = "TI_3D\\fold_convoluted52_SMALL.sgems";

	//string tiPath = "TI_3D\\checker_TI.SGEMS";
	string tiPath = "TI_3D\\fold_categorical.SGEMS";
	string seismicPath = "TI_3D\\fold_convoluted32.sgems";
	string syntheticPath = "TI_3D\\fold_convoluted52.sgems";

	short*** trainingImage;
	double*** seismicConditioning;
	double*** syntheticData;


	double thresholdMaxDistance = 0.5;
	int maxAmountOfPatternsInTheQueue = 100;
	int maxCandidates = 10;


	string realizationFileName1 = "Realizations\\myRealization.sgems";

	

};

#endif // !CONFIG_FILE_READER