#include <SDKDDKVer.h>
#include <stdio.h>
#include <tchar.h>
#include <unordered_map>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "ConfigFile.h"
#include "ToolBox.h"
#include "Coord.h"
#include "CPUTime.h"
#include "KMeansClustering.h"
#include "KMeans3DClustering.h"

using namespace std;

//void PrintShit(ConfigFile* configFile, short*** trainingImage)
//{
//
//	string arqName = "Realizations3D\\";
//
//	arqName = arqName + "Simulation_" + to_string(0) + "_" //+ configFile->GetRealizationFileName()
//		+ "R" + to_string(configFile->GetRealizationWidth()) + "x" + to_string(configFile->GetRealizationHeight()) + "x" + to_string(configFile->GetRealizationDepth()) + "_"
//		+ "T" + to_string(configFile->GetTemplateWidth()) + "x" + to_string(configFile->GetTemplateHeight()) + "x" + to_string(configFile->GetTemplateDepth()) + "_"
//		+ "O" + to_string(configFile->GetOverlapWidth()) + "x" + to_string(configFile->GetOverlapHeight()) + "x" + to_string(configFile->GetOverlapDepth()) //   +"_"
//		+ ".sgems";
//
//
//	//FILE* realizationFile = fopen(arqName.c_str(), "w");
//	ofstream myfile;
//	myfile.open(arqName);
//
//	//fprintf(realizationFile, "%d %d %d\n", configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
//	myfile << configFile->GetRealizationWidth() << " " << configFile->GetRealizationHeight() << " " << configFile->GetRealizationDepth() << endl;
//	//fprintf(realizationFile, "1\nfacies\n");
//	myfile << "1" << endl << "nfacies" << endl;
//
//	for (int d = 0; d < configFile->GetTIDepth(); d++)
//		for (int h = 0; h < configFile->GetTIHeight(); h++)
//			for (int w = 0; w < configFile->GetTIWidth(); w++)
//			{
//				//fprintf(realizationFile, "%d\n", simulationGrid[w][h][d]);
//				myfile << trainingImage[w][h][d] << endl;
//			}
//
//	//fclose(realizationFile);
//	myfile.close();
//
//
//
//}
//void PrintShit(ConfigFile* configFile, double*** trainingImage)
//{
//
//	string arqName = "Realizations3D\\";
//
//	arqName = arqName + "Simulation_" + to_string(0) + "_" //+ configFile->GetRealizationFileName()
//		+ "R" + to_string(configFile->GetRealizationWidth()) + "x" + to_string(configFile->GetRealizationHeight()) + "x" + to_string(configFile->GetRealizationDepth()) + "_"
//		+ "T" + to_string(configFile->GetTemplateWidth()) + "x" + to_string(configFile->GetTemplateHeight()) + "x" + to_string(configFile->GetTemplateDepth()) + "_"
//		+ "O" + to_string(configFile->GetOverlapWidth()) + "x" + to_string(configFile->GetOverlapHeight()) + "x" + to_string(configFile->GetOverlapDepth()) //   +"_"
//		+ ".sgems";
//
//
//	//FILE* realizationFile = fopen(arqName.c_str(), "w");
//	ofstream myfile;
//	myfile.open(arqName);
//
//	//fprintf(realizationFile, "%d %d %d\n", configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
//	myfile << configFile->GetTIWidth() << " " << configFile->GetTIHeight() << " " << configFile->GetTIDepth() << endl;
//	//fprintf(realizationFile, "1\nfacies\n");
//	myfile << "1" << endl << "nfacies" << endl;
//
//	for (int d = 0; d < configFile->GetTIDepth(); d++)
//		for (int h = 0; h < configFile->GetTIHeight(); h++)
//			for (int w = 0; w < configFile->GetTIWidth(); w++)
//			{
//				//fprintf(realizationFile, "%d\n", simulationGrid[w][h][d]);
//				myfile << trainingImage[w][h][d] << endl;
//			}
//
//	//fclose(realizationFile);
//	myfile.close();
//
//
//
//}
	//PrintShit(configFile, trainingImage);
	//PrintShit(configFile, syntheticSeismicData);
	//PrintShit(configFile, seismicConditioning);

//void ReduceSizeFromTIAndOthers()
//{
//	char* configPath = "Config.txt";
//	ConfigFile* configFile = new ConfigFile(configPath);
//
//	int originalWidth = 0;
//	int originalHeight = 0;
//	int originalDepth = 0;
//	string originalPath = "TI_3D\\fold_categorical.sgems";
//
//	int targetWidth = 52;
//	int targetHeight = 52;
//	int targetDepth = 52;
//
//	double*** originalImage = UniversalSGEMSreader::ReadSGEMS_3D_double_FAST(originalPath, originalWidth, originalHeight, originalDepth);
//
//	double*** newSizedMatrix = ToolBox::InitializeMatrix_double(targetWidth, targetHeight, targetDepth);
//
//	ConfigFile::GetSubMatrix(originalImage, 0, 0, 0, newSizedMatrix, targetWidth, targetHeight, targetDepth);
//
//	ToolBox::SaveSimulationToSGEMS(newSizedMatrix, configFile, 0);
//
//}


void RunWithoutMeasuringTime()
{
	char* configPath = "Config.txt";
	ConfigFile* configFile = new ConfigFile(configPath);								// load the Configuration file
	cout << "Loaded Configuration File" << endl;

	short*** trainingImage = configFile->GetTrainingImage();							// load the Training Image (uses singleton methodology)
	cout << "Loaded Training Image" << endl;

	double*** syntheticSeismicData = configFile->GetSyntheticSeismicData();				// load the Synthetic Seismic Data (uses Singleton methodology)
	cout << "Loaded Seismic Data" << endl;

	double*** seismicConditioning = configFile->GetSeismicConditioning();				// load the Seismic Conditioning Data (uses singleton methodology)
	cout << "Loaded Seismic Conditioning" << endl;

	short*** simulationGrid = ToolBox::InitializeMatrix_short(configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	ToolBox::ResetMatrix(simulationGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());

	vector<Coord> simulationPath	= *ToolBox::DefineSimulationPath(configFile);		// this guy will hold the Simulation Path
	cout << "Created Simulation Path" << endl;

	vector<Coord> patternsPath		= *ToolBox::DefinePatternsPath(configFile);			// this guy will hold the coordinates of all existing Patterns.
	cout << "Created Patterns Path" << endl;


	// Here I give some margin to find which type of Clustering should be done over the Patterns from the Training Image, gotten from the configFile->
	vector<Coord> optimizedPatternCoords = *ToolBox::OptimizePatternCoordinates(trainingImage, patternsPath, configFile);
	cout << "Optimized the Patterns" << endl;


	// Here I give some margin to find which type of Clustering should be done over the Seismic Data Events (SeismCond + SimPath), gotten from the configFile->
	vector<Coord> optimizedSeismicDataEventsCoords = *ToolBox::OptimizeSeismicCoordinates(seismicConditioning, simulationPath, configFile);
	cout << "Optimized the Seismic Data" << endl;


	
	// Associate the Patterns from the TI to the Seismic Data Events, by comparisons between the SDE and the Synthetic Seismic Data that corresponds to the PTI.
	unordered_map<Coord, vector<Coord>*, Coord> association = *ToolBox::Associate(
		optimizedSeismicDataEventsCoords, seismicConditioning,			// the Optimized Seismic Data Events,
		optimizedPatternCoords, syntheticSeismicData,					// the Synthetic Seismic Data (represented by Optimized Patterns from the Training Image),
		configFile);													// and the support data.
	cout << "Association between Patterns and Seismic Data Events created" << endl;



	// Now comes the easyest part: the main loop of simulations!!!
	int amountOfSimulations = configFile->GetAmountOfSimulations();
	for (int currentSimulation = 0; currentSimulation < amountOfSimulations; currentSimulation++)
	{
		cout << "		STARTING SIMULATION " << (currentSimulation+1) << endl << endl;

		ToolBox::RunSimulationWithSeismicConditioning(simulationGrid, simulationPath, association, seismicConditioning, trainingImage, configFile);
		cout << "	Simulation finished!" << endl;							// so here is where the	WHOLE simulation process takes place...

		ToolBox::SaveSimulationToSGEMS(simulationGrid, configFile, currentSimulation);
		cout << "	Simulation saved!" << endl;								// here we have saved it to a new file...

		ToolBox::ResetMatrix(simulationGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
		cout << "	Simulation grid was reset!" << endl;					// and here, we reset the Simulation Grid. 

		cout << endl << endl;
	}
	
	ToolBox::DeleteMatrix(simulationGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
}


void RunMeasuringTime(int repetitions)
{
	char* configPath = "Config.txt";
	ConfigFile* configFile = new ConfigFile(configPath);								// load the Configuration file
	cout << "Loaded Configuration File" << endl;

	short*** trainingImage = configFile->GetTrainingImage();							// load the Training Image (uses singleton methodology)
	cout << "Loaded Training Image" << endl;

	double*** syntheticSeismicData = configFile->GetSyntheticSeismicData();				// load the Synthetic Seismic Data (uses Singleton methodology)
	cout << "Loaded Seismic Data" << endl;

	double*** seismicConditioning = configFile->GetSeismicConditioning();				// load the Seismic Conditioning Data (uses singleton methodology)
	cout << "Loaded Seismic Conditioning" << endl;

	short*** simulationGrid = ToolBox::InitializeMatrix_short(configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	ToolBox::ResetMatrix(simulationGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());


	double time1 = 0;


	cout << "Here comes the repetitions: we will perform " << repetitions << " repetitions." << endl << endl;
	for (int currentIteration = 0; currentIteration < repetitions; currentIteration++)
	{
		cout << "	Iteration #" << (currentIteration + 1) << endl;

		time1 = get_cpu_time();

		vector<Coord> simulationPath = *ToolBox::DefineSimulationPath(configFile);		// this guy will hold the Simulation Path
		cout << "Created Simulation Path" << endl;

		vector<Coord> patternsPath = *ToolBox::DefinePatternsPath(configFile);			// this guy will hold the coordinates of all existing Patterns.
		cout << "Created Patterns Path" << endl;


		// Here I give some margin to find which type of Clustering should be done over the Patterns from the Training Image, gotten from the configFile->
		vector<Coord> optimizedPatternCoords = *ToolBox::OptimizePatternCoordinates(trainingImage, patternsPath, configFile);
		cout << "Optimized the Patterns" << endl;


		// Here I give some margin to find which type of Clustering should be done over the Seismic Data Events (SeismCond + SimPath), gotten from the configFile->
		vector<Coord> optimizedSeismicDataEventsCoords = *ToolBox::OptimizeSeismicCoordinates(seismicConditioning, simulationPath, configFile);
		cout << "Optimized the Seismic Data" << endl;



		// Associate the Patterns from the TI to the Seismic Data Events, by comparisons between the SDE and the Synthetic Seismic Data that corresponds to the PTI.
		unordered_map<Coord, vector<Coord>*, Coord> association = *ToolBox::Associate(
			optimizedSeismicDataEventsCoords, seismicConditioning,			// the Optimized Seismic Data Events,
			optimizedPatternCoords, syntheticSeismicData,					// the Synthetic Seismic Data (represented by Optimized Patterns from the Training Image),
			configFile);													// and the support data.
		cout << "Association between Patterns and Seismic Data Events created" << endl;



		// Now comes the easyest part: the main loop of simulations!!!
		int amountOfSimulations = configFile->GetAmountOfSimulations();
		for (int currentSimulation = 0; currentSimulation < amountOfSimulations; currentSimulation++)
		{
			cout << "		STARTING SIMULATION " << (currentSimulation + 1) << endl << endl;

			ToolBox::RunSimulationWithSeismicConditioning(simulationGrid, simulationPath, association, seismicConditioning, trainingImage, configFile);
			cout << "	Simulation finished!" << endl;							// so here is where the	WHOLE simulation process takes place...

			ToolBox::SaveSimulationToSGEMS(simulationGrid, configFile, currentSimulation);
			cout << "	Simulation saved!" << endl;								// here we have saved it to a new file...

			ToolBox::ResetMatrix(simulationGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
			cout << "	Simulation grid was reset!" << endl;					// and here, we reset the Simulation Grid. 

			cout << endl << endl;
		}


	}



	ToolBox::DeleteMatrix(simulationGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
}

void readData2D(double **matrix) {
	
	ifstream myfile("data.txt");
	if (myfile.is_open())
	{
		int j = 0;
		string line;
		while (getline(myfile, line))
		{
			double a, b, c, d;
			cout << "checking line " << j + 1 << ": ";
			istringstream iss(line);
			iss >> a >> b >> c >> d;
			matrix[0][j] = a;
			matrix[1][j] = b;
			matrix[2][j] = c;
			matrix[3][j] = d;
			cout << matrix[0][j] << " " << matrix[1][j] << " " << matrix[2][j]  << " " << matrix[3][j] << endl;
			
			j++;
		}
		cout << "file read" << endl;
		myfile.close();
	}

	else cout << "Unable to open file";
		
}

void readData3D(double ***matrix) {

	ifstream myfile("data.txt");
	if (myfile.is_open())
	{
		int j = 0;
		string line;
		while (getline(myfile, line))
		{
			double a, b, c, d;
			cout << "checking line " << j + 1 << ": ";
			istringstream iss(line);
			iss >> a >> b >> c >> d;
			matrix[0][0][j] = a;
			matrix[1][0][j] = b;
			matrix[0][1][j] = c;
			matrix[1][1][j] = d;
			cout << matrix[0][0][j] << " " << matrix[1][0][j] << " " << matrix[0][1][j] << " " << matrix[1][1][j] << endl;

			j++;
		}
		cout << "file read" << endl;
		myfile.close();
	}

	else cout << "Unable to open file";

}

void RunKMeansClusteringTest() {
	int K = 3;
	int Max_interations = 50;
	cout << "-----------starting KMEANS 2D test------------" << endl;
	ToolBox* toolbox = new ToolBox();
	KMeansClustering<double>* kmeans = new KMeansClustering<double>(K, Max_interations);
	double **matrix = toolbox->InitializeMatrix_double(4, 150);
	readData2D(matrix);
	kmeans->run(matrix, 4, 150, 4, 1);

}

void RunKMeans3DClusteringTest() {
	int K = 3;
	int Max_interations = 50;
	cout << "-----------starting KMEANS 3D test------------" << endl;
	ToolBox* toolbox = new ToolBox();
	KMeansClustering3D<double>* kmeans = new KMeansClustering3D<double>(K, Max_interations);
	double ***matrix = toolbox->InitializeMatrix_double(2, 2, 150);

	readData3D(matrix);
	kmeans->run(matrix, 2, 2, 150, 2, 2, 1);


}

int main()
{
	RunKMeans3DClusteringTest();
	//RunKMeansClusteringTest();
	//RunWithoutMeasuringTime();
	//RunMeasuringTime(10);

    return 0;																// finally, cleanup and return bro.
}