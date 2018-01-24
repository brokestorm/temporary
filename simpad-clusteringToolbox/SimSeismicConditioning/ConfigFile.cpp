#include "ConfigFile.h"


ConfigFile::ConfigFile(char* path)
{
	//this->path = path;


	ifstream fs(path);

	string line = "";
	vector<string> tokens;

	line = FileReadingToolBox::GetNextValidLine(fs);

	//for (string line; getline(fs, line); ){}
	do
	{

		// process the line....

		//vector<string> tokens = ToolBox::SplitString(line, ":");




		line = FileReadingToolBox::GetNextValidLine(fs);
	} while (!line.empty());




	//fscanf(configFile, "Training Image File: %s\n", TIFile);
	//fscanf(configFile, "Training Image Size: %dx%dx%d\n", &TIWidth, &TIHeight, &TIDepth);

	//fscanf(configFile, "Realization File Name: %s\n", realizationFile);
	//fscanf(configFile, "Realization Size: %dx%dx%d\n", &realizationWidth, &realizationHeight, &realizationDepth);

	//fscanf(configFile, "Template Size : %dx%dx%d\n", &initialTemplateWidth, &initialTemplateHeight, &initialTemplateDepth);
	//fscanf(configFile, "Overlap: %dx%dx%d\n", &overlapWidth, &overlapHeight, &overlapDepth);

	//fscanf(configFile, "Simulation Path: %s\n", simulationPath);
	//fscanf(configFile, "Similarity Function: %s\n", similarityFunction);
	//fscanf(configFile, "Random seed: %d\n", &randomSeed);
	//fscanf(configFile, "Amount of Simulations: %d\n", &numOfSims);

	//fscanf(configFile, "Training Image Convoluted File : %s\n", TI_CONV_FILE);			//TI_3D\fold_seismic32.sgems
	//fscanf(configFile, "Seismic Conditioning File : %s\n", SEISM_FILE);									// TI_3D\fold_seismic32.sgems

	//FILE* configFile;
	//configFile = fopen(path, "r");			// "Config.txt"
	//if (configFile == NULL)
	//{
	//	printf("please choose a valid config file\n");
	//	printf("press enter to abort\n");
	//	getchar();
	//	exit(0);
	//}
	//fscanf(configFile, "Number of Facies : %d\n", &numFacies);
	//
	//
	//
	//
	//
	//fscanf(configFile, "Hard Data File: %s\n", HDFile);
	////fscanf(configFile, "Final Template Size : %dx%dx%d\n", &maxTemplateWidth, &maxTemplateHeight, &maxTemplateDepth);
	////fscanf(configFile, "Template Size Increment : %dx%dx%d\n", &templateWidthIncrement, &templateHeightIncrement, &templateDepthIncrement);
	//fscanf(configFile, "Compression: %s\n", compression);
	//fscanf(configFile, "ANN L: %d\n", &ANNL);
	//fscanf(configFile, "ANN K: %d\n", &ANNK);
	//fscanf(configFile, "Alpha: %f\n", &alpha);
	//fscanf(configFile, "Max Candidates: %d\n", &maxCandidates);
	//fscanf(configFile, "Weight Matrix : %s\n", WMFile);

	//if (randomSeed == 0)
	//{
	//	randomSeed = time(NULL);
	//}
}

ConfigFile::~ConfigFile()
{
}



int ConfigFile::GetTIWidth()
{
	return tiWidth;
}
int ConfigFile::GetTIHeight()
{
	return tiHeight;
}
int ConfigFile::GetTIDepth()
{
	return tiDepth;
}

int ConfigFile::GetTemplateWidth()
{
	return templateWidth;
}
int ConfigFile::GetTemplateHeight()
{
	return templateHeight;
}
int ConfigFile::GetTemplateDepth()
{
	return templateDepth;
}

int ConfigFile::GetRealizationWidth()
{
	return realizationWidth;
}
int ConfigFile::GetRealizationHeight()
{
	return realizationHeight;
}
int ConfigFile::GetRealizationDepth()
{
	return realizationDepth;
}

int ConfigFile::GetOverlapWidth()
{
	return overlapWidth;
}
int ConfigFile::GetOverlapHeight()
{
	return overlapHeight;
}
int ConfigFile::GetOverlapDepth()
{
	return overlapDepth;
}

int ConfigFile::GetAmountOfSimulations()
{
	return amountOfSimulations;
}
double ConfigFile::GetMinimumCompatibility()
{
	return minimumCompatibility;
}
double ConfigFile::GetSigmaSquared()
{
	return sigmaSquared;
}



short*** ConfigFile::GetTrainingImage()
{
	// SINGLETON IMPLEMENTATION!!!
	if (trainingImage == nullptr)
	{
		int auxW = 0;
		int auxH = 0;
		int auxD = 0;
		auto originalImage = UniversalSGEMSreader::ReadSGEMS_3D_short_FAST(tiPath, auxW, auxH, auxD);		// here the file was read and the ACTUAL sizes initialized.
		trainingImage = InitializeMatrix_short(tiWidth, tiHeight, tiDepth);
		GetSubMatrix(originalImage, 0, 0, 0, trainingImage, tiWidth, tiHeight, tiDepth);
		DeleteMatrix(originalImage, auxW, auxH, auxD);
		for (int d = 0; d < tiDepth; d++)
			for (int h = 0; h < tiHeight; h++)
				for (int w = 0; w < tiWidth; w++)
					trainingImage[w][h][d]++;

	}
		return trainingImage;
}
double*** ConfigFile::GetSyntheticSeismicData()
{
	// SINGLETON IMPLEMENTATION!!!
	if (syntheticData == nullptr)
	{
		int auxW = 0;
		int auxH = 0;
		int auxD = 0;
		auto aux = UniversalSGEMSreader::ReadSGEMS_3D_double_FAST(syntheticPath, auxW, auxH, auxD);
		
		syntheticData = InitializeMatrix_double(tiWidth, tiHeight, tiDepth);
		GetSubMatrix(aux, 0, 0, 0, syntheticData, tiWidth, tiHeight, tiDepth);
		DeleteMatrix(aux, auxW, auxH, auxD);
	}
	return syntheticData;
}
double*** ConfigFile::GetSeismicConditioning()
{
	// SINGLETON IMPLEMENTATION!!!
	if (seismicConditioning == nullptr)
	{
		int auxW = 0;
		int auxH = 0;
		int auxD = 0;
		auto aux = UniversalSGEMSreader::ReadSGEMS_3D_double_FAST(seismicPath, auxW, auxH, auxD);

		seismicConditioning = InitializeMatrix_double(realizationWidth, realizationHeight, realizationDepth);
		GetSubMatrix(aux, 0, 0, 0, seismicConditioning, realizationWidth, realizationHeight, realizationDepth);
		DeleteMatrix(aux, auxW, auxH, auxD);
	}
	return seismicConditioning;
}


double ConfigFile::GetThresholdMaxDistance()
{
	return thresholdMaxDistance;
}
int ConfigFile::GetMaxAmountOfPatternsInTheQueue()
{
	return maxAmountOfPatternsInTheQueue;
}
int ConfigFile::GetMaxCandidates()
{
	return maxCandidates;
}



string ConfigFile::GetRealizationFileName11()
{
	return realizationFileName1;
}


template <class T>
void ConfigFile::GetSubMatrix(T*** bigMatrix, int x, int y, int z, T*** smallMatrix, int matrixWidth, int matrixHeight, int matrixDepth)
{
	for (int d = 0; d < matrixDepth; d++)
		for (int h = 0; h < matrixHeight; h++)
			for (int w = 0; w < matrixWidth; w++)
				smallMatrix[w][h][d] = bigMatrix[x + w][y + h][z + d];

}

short*** ConfigFile::InitializeMatrix_short(int width, int height, int depth)
{
	short*** matrix = new short**[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new short*[height];
		for (int j = 0; j < height; j++)
		{
			matrix[i][j] = new short[depth];
		}
	}

	return matrix;
}

double*** ConfigFile::InitializeMatrix_double(int width, int height, int depth)
{
	double*** matrix = new double**[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new double*[height];
		for (int j = 0; j < height; j++)
		{
			matrix[i][j] = new double[depth];
		}
	}

	return matrix;
}

void ConfigFile::DeleteMatrix(short*** matrix, int width, int height, int depth)
{
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			//for (int k = 0; k < depth; k++)
			//{
			//	delete matrix[i][j][k];
			//}
			delete[] matrix[i][j];

		}
		delete[] matrix[i];
	}

	delete[] matrix;
}

void ConfigFile::DeleteMatrix(double*** matrix, int width, int height, int depth)
{
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			//for (int k = 0; k < depth; k++)
			//{
			//	delete matrix[i][j][k];
			//}
			delete[] matrix[i][j];

		}
		delete[] matrix[i];
	}

	delete[] matrix;
}

