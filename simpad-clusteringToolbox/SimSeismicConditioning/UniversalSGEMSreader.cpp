#include "UniversalSGEMSreader.h"

short*** UniversalSGEMSreader::ReadSGEMS_3D_short1(string path, int& auxW, int& auxH, int& auxD)
{
	ifstream fs(path);					// initialize the IFSTREAM to read the file.

	string line = FileReadingToolBox::GetNextValidLine(fs);						// the dimensions, as a raw string
	vector<string> splittedLine = FileReadingToolBox::SplitString(line, " ");	// the dimensions, as separated strings.

	if (splittedLine.size() == 3)
	{
		int tiWidth = stoi(splittedLine[0]);										// parsed the Width
		int tiHeight = stoi(splittedLine[1]);										// parsed the Height
		int tiDepth = stoi(splittedLine[2]);										// parsed the Depth

		short*** matrix = InitializeMatrix_short(tiWidth, tiHeight, tiDepth);	// now, here, initialize the Matrix.


		int amountOfCols = stoi(FileReadingToolBox::GetNextValidLine(fs));			// get and parse the amount of columns (should be only one)
		vector<string> namesOfCols;													// for just in case, let's store them in here
		for (int currentCol = 0; currentCol < amountOfCols; currentCol++)
			namesOfCols.push_back(FileReadingToolBox::GetNextValidLine(fs));		// yep, storing them in the vector.



		for (int d = 0; d < tiDepth; d++)
		{
			for (int h = 0; h < tiHeight; h++)
			{
				for (int w = 0; w < tiWidth; w++)
				{
					line = FileReadingToolBox::GetNextValidLine(fs);				// get each line of the data...
					splittedLine = FileReadingToolBox::SplitString(line, " ");		// split the data (in case it should be split)

					matrix[w][h][d] = (short)stoi(splittedLine[0]);					// PROBLEM!!! IF THERE'S MORE THAN 1 COLUMN, I LOSE THE REST!!!
				}
			}
		}

		auxW = tiWidth;
		auxH = tiHeight;
		auxD = tiDepth;

		return matrix;
	}
	else if(splittedLine.size() == 1)
	{
		int amountOfCols = stoi(FileReadingToolBox::GetNextValidLine(fs));			// get and parse the amount of columns (should be only one)
		vector<string> namesOfCols;													// for just in case, let's store them in here
		for (int currentCol = 0; currentCol < amountOfCols; currentCol++)
			namesOfCols.push_back(FileReadingToolBox::GetNextValidLine(fs));		// yep, storing them in the vector.
		
	}
	else
	{
		return nullptr;
	}
}

short*** UniversalSGEMSreader::ReadSGEMS_3D_short_FAST(string path, int& auxW, int& auxH, int& auxD)
{
	ifstream fs(path);					// initialize the IFSTREAM to read the file.

	string line = FileReadingToolBox::GetNextValidLine(fs);						// the dimensions, as a raw string
	vector<string> splittedLine = FileReadingToolBox::SplitString(line, " ");	// the dimensions, as separated strings.
	

	int tiWidth = stoi(splittedLine[0]);										// parsed the Width
	int tiHeight = stoi(splittedLine[1]);										// parsed the Height
	int tiDepth = stoi(splittedLine[2]);										// parsed the Depth

	short*** matrix = InitializeMatrix_short(tiWidth, tiHeight, tiDepth);	// now, here, initialize the Matrix.


	int amountOfCols = stoi(FileReadingToolBox::GetNextValidLine(fs));			// get and parse the amount of columns (should be only one)
	vector<string> namesOfCols;													// for just in case, let's store them in here
	for (int currentCol = 0; currentCol < amountOfCols; currentCol++)
		namesOfCols.push_back(FileReadingToolBox::GetNextValidLine(fs));		// yep, storing them in the vector.



	for (int d = 0; d < tiDepth; d++)
	{
		for (int h = 0; h < tiHeight; h++)
		{
			for (int w = 0; w < tiWidth; w++)
			{
				line = FileReadingToolBox::GetNextValidLine(fs);				// get each line of the data...
				//splittedLine = FileReadingToolBox::SplitString(line, " ");		// split the data (in case it should be split)

				matrix[w][h][d] = (short)stoi(line);					// PROBLEM!!! IF THERE'S MORE THAN 1 COLUMN, I LOSE THE REST!!!
			}
		}
	}

	auxW = tiWidth;
	auxH = tiHeight;
	auxD = tiDepth;

	return matrix;
	
}


double*** UniversalSGEMSreader::ReadSGEMS_3D_double_FAST(string path, int& auxW, int& auxH, int& auxD)
{
	ifstream fs(path);					// initialize the IFSTREAM to read the file.

	string line = FileReadingToolBox::GetNextValidLine(fs);						// the dimensions, as a raw string
	vector<string> splittedLine = FileReadingToolBox::SplitString(line, " ");	// the dimensions, as separated strings.


	int tiWidth = stoi(splittedLine[0]);										// parsed the Width
	int tiHeight = stoi(splittedLine[1]);										// parsed the Height
	int tiDepth = stoi(splittedLine[2]);										// parsed the Depth
	double*** matrix = InitializeMatrix_double(tiWidth, tiHeight, tiDepth);	// now, here, initialize the Matrix.


	int amountOfCols = stoi(FileReadingToolBox::GetNextValidLine(fs));			// get and parse the amount of columns (should be only one)
	vector<string> namesOfCols;													// for just in case, let's store them in here
	for (int currentCol = 0; currentCol < amountOfCols; currentCol++)
		namesOfCols.push_back(FileReadingToolBox::GetNextValidLine(fs));		// yep, storing them in the vector.



	for (int d = 0; d < tiDepth; d++)
	{
		for (int h = 0; h < tiHeight; h++)
		{
			for (int w = 0; w < tiWidth; w++)
			{
				line = FileReadingToolBox::GetNextValidLine(fs);				// get each line of the data...
				//splittedLine = FileReadingToolBox::SplitString(line, " ");		// split the data (in case it should be split)

				matrix[w][h][d] = stod(line);					// PROBLEM!!! IF THERE'S MORE THAN 1 COLUMN, I LOSE THE REST!!!
			}
		}
	}

	
	auxW = tiWidth;
	auxH = tiHeight;
	auxD = tiDepth;
	return matrix;
}


short*** UniversalSGEMSreader::InitializeMatrix_short(int width, int height, int depth)
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


double*** UniversalSGEMSreader::InitializeMatrix_double(int width, int height, int depth)
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