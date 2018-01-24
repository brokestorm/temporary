#ifndef _TOOLBOX
#define _TOOLBOX

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <windows.h>
#include <queue>

#include "ConfigFile.h"
#include "Coord.h"


using namespace std;

template <class T>
struct LessThanByFirst
{
	bool operator()(const pair<double, T>& lhs, const pair<double, T>& rhs) const
	{
		return lhs.first < rhs.first;
	}
};

template <class T>
struct GreaterThanByFirst
{
	bool operator()(const pair<double, T>& lhs, const pair<double, T>& rhs) const
	{
		return lhs.first > rhs.first;
	}
};

//struct LessThanByFirst
//{
//	bool operator()(const pair<double, int>& lhs, const pair<double, int>& rhs) const
//	{
//		return lhs.first < rhs.first;
//	}
//};
//struct GreaterThanByFirst
//{
//	bool operator()(const pair<double, int>& lhs, const pair<double, int>& rhs) const
//	{
//		return lhs.first > rhs.first;
//	}
//};


class ToolBox
{
public:

	static vector<Coord>* DefineSimulationPath(ConfigFile* configFile);
	static vector<Coord>* DefinePatternsPath(ConfigFile* configFile);

	static vector<Coord>* OptimizeSeismicCoordinates(double*** seismicConditioning, vector<Coord>& simulationPath, ConfigFile* configFile);
	static vector<Coord>* OptimizePatternCoordinates(short*** trainingImage, vector<Coord>& patternsPath, ConfigFile* configFile);

	// Associate the Patterns from the TI to the Seismic Data Events, by comparisons between the SDE and the Synthetic Seismic Data that corresponds to the PTI.
	static unordered_map<Coord, vector<Coord>*, Coord>* Associate(
		vector<Coord>& optimizedSeismicDataEventsCoords,
		double*** seismicConditioning,
		vector<Coord>& optimizedPatternsCoords,
		double*** syntheticSeismicData,
		ConfigFile* configFile);		// and the support data.

	static void ToolBox::RunSimulationWithSeismicConditioning(
		short*** simulationGrid,
		vector<Coord>& simulationPath,
		unordered_map<Coord, vector<Coord>*, Coord>& association,
		double*** seismicConditioning,
		short*** trainingImage,
		ConfigFile* configFile);

	static void ToolBox::RunSimulationWithSeismicConditioning_SeismicImage(
		short*** simulationGrid,
		vector<Coord>& simulationPath,
		unordered_map<Coord, vector<Coord>*, Coord>& association,
		double*** seismicConditioning,
		short*** trainingImage,
		double*** syntheticSeismicData,
		int currentSimulation,
		ConfigFile* configFile);

	static short*** InitializeMatrix_short(int width, int height, int depth);
	static short*** InitializeMatrix_short(Coord& dimensions);

	static short**  InitializeMatrix_short(int width, int height);

	static int**  InitializeMatrix_int(int width, int height);

	static short*** InitializeMatrix_short_Random(int width, int height, int depth, int minRandom, int maxRandom);
	static short**  InitializeMatrix_short_Random(int width, int height, int minRandom, int maxRandom);

	static double*** InitializeMatrix_double(int width, int height, int depth);
	static double**  InitializeMatrix_double(int width, int height);

	static double*** InitializeMatrix_double_Random(int width, int height, int depth, int minRandom, int maxRandom);
	static double**  InitializeMatrix_double_Random(int width, int height, int minRandom, int maxRandom);


	//template <class T>
	//static void DeleteMatrix(T*** simulationGrid, int width, int height, int depth);
	static void DeleteMatrix(short*** matrix, Coord& dimensions);
	static void DeleteMatrix(short***  simulationGrid, int width, int height, int depth);
	static void DeleteMatrix(double*** simulationGrid, int width, int height, int depth);
	static void DeleteMatrix(double**  simulationGrid, int width, int height);
	static void DeleteMatrix(int**  simulationGrid, int width, int height);

	static void SaveSimulationToSGEMS(short*** simulationGrid, ConfigFile* configFile, int currentSimulation);
	static void SaveSimulationToSGEMS(double*** simulationGrid, ConfigFile* configFile, int currentSimulation);

	template <class T>
	static void ResetMatrix(T*** matrix, int width, int height, int depth);
	static void ResetMatrix(short*** matrix, int width, int height, int depth);
	//static void ResetMatrix(double*** matrix, int width, int height, int depth);

	static pair<int, int>* Depair_CantorPairing(int z);
	static int Pair_CantorPairing(int x, int y);
	static void TestPairingDepairing(int maxLimit);
	static bool Create_Directory(char* folderName);
	//static void GetDataEvent(double*** seismicConditioning, Coord seismicDataEvent, double*** seismicMatrix, ConfigFile* configFile);
	//static void GetDataEvent(short*** seismicConditioning, Coord seismicDataEvent, short*** seismicMatrix, ConfigFile* configFile);


	template <class T>
	static void GetSubMatrix(T*** bigMatrix, Coord& sizeOfBigMatrix, Coord& startInBigMatrix, T*** smallMatrix, Coord& sizeOfSmallMatrix);
	template <class T>
	static void PasteSubMatrix(T*** bigMatrix, Coord& sizeOfBigMatrix, Coord& startInBigMatrix, T*** smallMatrix, Coord& sizeOfSmallMatrix);
	//template <class T>
	//static void PasteSubMatrix_MEBC(T*** seismicConditioning, Coord seismicDataEvent, T*** seismicMatrix, ConfigFile* configFile);
	static void PasteSubMatrix_MEBC(short*** bigMatrix, Coord& sizeOfBigMatrix, Coord& startInBigMatrix, short*** smallMatrix, Coord& sizeOfSmallMatrix, ConfigFile* configFile);
	template <class T>
	static void PasteSubMatrix_AveragedOverlap(T*** seismicConditioning, Coord seismicDataEvent, T*** seismicMatrix, ConfigFile* configFile);

	static void FLIP_SIDEMATRIX_INTO_TOPMATRIX(short*** sideMatrix, Coord& sideMatrixDimensions, short*** topMatrix, Coord& topMatrixDimensions);
	static void FLIP_TOPMATRIX_INTO_SIDEMATRIX(short*** topMatrix, Coord& topMatrixDimensions, short*** sideMatrix, Coord& sideMatrixDimensions);
	static void FLIP_FRONTMATRIX_INTO_TOPMATRIX(short*** frontMatrix, Coord& frontMatrixDimensions, short*** topMatrix, Coord& topMatrixDimensions);
	static void FLIP_TOPMATRIX_INTO_FRONTMATRIX(short*** topMatrix, Coord& topMatrixDimensions, short*** frontMatrix, Coord& frontMatrixDimensions);



	static void SelectPattern(short*** trainingImage, vector<Coord>& associatedPatterns, short*** dataEventMatrix, short*** selectedPattern, ConfigFile* configFile);
	static Coord SelectPattern(short*** trainingImage, vector<Coord>& associatedPatterns, short*** dataEventMatrix, ConfigFile* configFile);


	static double CalculateCompatibility(double*** seismicMatrix, double*** syntheticMatrix, int templateWidth, int templateHeight, int templateDepth, double sigmaSquared);
	static double EuclideanDistance(double*** matrixA, double*** matrixB, int width, int height, int depth);
	static double EuclideanDistance(short*** matrixA, short*** matrixB, int width, int height, int depth);
	static double EuclideanDistance_WithoutZeroes(short*** matrixA, short*** matrixB, int width, int height, int depth);

	//template <class T>
	static void PrintMatrix(short** matrix, int width, int height);
	static void PrintMatrix(double** matrix, int width, int height);




	static int IndexOfSmallestValue(double* vector, int vectorSize);

	template<class T>
	static void InitializeSurfaceError_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, double** surfaceError);
	template<class T>
	static void InitializeSurfaceError_3D(T*** AMatrix, T*** BMatrix, int sizeX_template, int sizeY_template, int sizeZ_overlap, double*** surfaceError);

	static void InitializeCumulativeMinimumErrorSurface_2D(double** surfaceError, int sizeX_template, int sizeY_overlap, double** cumulMinError);
	static void InitializeCumulativeMinimumErrorSurface_3D(double*** surfaceError, int sizeX_template, int sizeY_template, int sizeZ_overlap, double*** cumulMinError);

	static void FindBoundaryHairline(double** cumulMinError, int sizeX_template, int sizeY_overlap, int* boundary);
	static void FindBoundarySurface(double*** cumulMinError, int sizeX_template, int sizeY_template, int sizeZ_overlap, int** boundary);

	template<class T>
	static void MergeMatricesUsingBoundaryCut_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, int* boundary, T** answer);
	template<class T>
	static void MergeMatricesUsingBoundaryCut_3D(T*** AMatrix, T*** BMatrix, int sizeX_template, int sizeY_template, int sizeZ_overlap, int** boundary, T*** answer);

	template<class T>
	static void MEBC_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, T** answer);

	template<class T>
	static void MEBC_3D(T*** AMatrix, T*** BMatrix, int sizeX_template, int sizeY_template, int sizeZ_overlap, T*** answer);


private:

	static unordered_map<Coord, vector<Coord>*, Coord>* ToolBox::AssociateByPercentageOfDistance(
		vector<Coord>& seismicDataEventsCoords,
		double*** seismicConditioning,
		vector<Coord>& patternCoords,
		double*** syntheticSeismicData,
		ConfigFile* configFile);

	static unordered_map<Coord, vector<Coord>*, Coord>* ToolBox::AssociateByPercentageOfTheTotal(
		vector<Coord>& seismicDataEventsCoords,
		double*** seismicConditioning,
		vector<Coord>& patternCoords,
		double*** syntheticSeismicData,
		ConfigFile* configFile);



};
#endif // !_TOOLBOX





//template <class T>
//static void MEBC_3D(T*** AMatrix, T*** BMatrix, int sizeX, int sizeY, int sizeZ, T*** answer);
//template <class T>
//static T** MEBC_2D(T** AMatrix, T** BMatrix, int sizeX, int sizeY, T** answer);// , ConfigFile* configFile);
//template <class T>
//static void PasteSubMatrix_MEBC_Mahmud(T*** bigMatrix, Coord startInBigMatrix, T*** smallMatrix, ConfigFile* configFile);
//template <class T>
//static void PasteSubMatrix_MEBC(T*** seismicConditioning, Coord seismicDataEvent, T*** seismicMatrix, ConfigFile* configFile);
