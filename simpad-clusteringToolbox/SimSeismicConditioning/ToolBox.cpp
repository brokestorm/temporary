#include "ToolBox.h"

#pragma region Definition of Simulation and Patterns Path

vector<Coord>* ToolBox::DefineSimulationPath(ConfigFile* configFile)
{
	vector<Coord>* simulationPath = new vector<Coord>();					// initializing some data here...

	int realizationWidth = configFile->GetRealizationWidth();				// initializing some data here...
	int realizationHeight = configFile->GetRealizationHeight();				// initializing some data here...
	int realizationDepth = configFile->GetRealizationDepth();				// initializing some data here...

	int templateWidth = configFile->GetTemplateWidth();					// initializing some data here...
	int templateHeight = configFile->GetTemplateHeight();					// initializing some data here...
	int templateDepth = configFile->GetTemplateDepth();					// initializing some data here...

	int overlapWidth = configFile->GetOverlapWidth();						// initializing some data here...
	int overlapHeight = configFile->GetOverlapHeight();						// initializing some data here...
	int overlapDepth = configFile->GetOverlapDepth();						// initializing some data here...

	int t_oW = templateWidth - overlapWidth;								// initializing some data here...
	int t_oH = templateHeight - overlapHeight;								// initializing some data here...
	int t_oD = templateDepth - overlapDepth;								// initializing some data here...

	int maxW_SG = realizationWidth / t_oW;									// + (realizationWidth  % t_oW == 0) ? 0 : 1;
	int maxH_SG = realizationHeight / t_oH;									// + (realizationHeight % t_oH == 0) ? 0 : 1;
	int maxD_SG = realizationDepth / t_oD;									// + (realizationDepth  % t_oD == 0) ? 0 : 1;

	Coord currentCoord;														// and this guy will hold the current coordinate.

	for (int k = 0; k < maxD_SG; k++)										// foreach Data Event coordinate in the Simulation Grid (and Seismic Conditioning)
		for (int j = 0; j < maxH_SG; j++)
			for (int i = 0; i < maxW_SG; i++)
				simulationPath->push_back(Coord(i*t_oW, j*t_oH, k*t_oD));	// get the current coordinate, and include it in the Simulation Path.

	return simulationPath;													// aaaaaaaaaaand finally, return the lovely path.
}

vector<Coord>* ToolBox::DefinePatternsPath(ConfigFile* configFile)
{
	vector<Coord>* patternsPath = new vector<Coord>();			// initializing some data here...

	int tiWidth = configFile->GetTIWidth();					// initializing some data here...
	int tiHeight = configFile->GetTIHeight();					// initializing some data here...
	int tiDepth = configFile->GetTIDepth();					// initializing some data here...

	int templateWidth = configFile->GetTemplateWidth();		// initializing some data here...
	int templateHeight = configFile->GetTemplateHeight();		// initializing some data here...
	int templateDepth = configFile->GetTemplateDepth();		// initializing some data here...

	int maxW_TI = tiWidth - templateWidth + 1;					// initializing some data here...
	int maxH_TI = tiHeight - templateHeight + 1;				// initializing some data here...
	int maxD_TI = tiDepth - templateDepth + 1;					// initializing some data here...

	Coord currentCoord;											// and this guy will hold the current coordinate.

	for (int k = 0; k < maxD_TI; k++)							// foreach Pattern coordinate in the Training Image...
		for (int j = 0; j < maxH_TI; j++)
			for (int i = 0; i < maxW_TI; i++)
				patternsPath->push_back(Coord(i, j, k));		// get the current coordinate, and include it in the Simulation Path.

	return patternsPath;										// aaaaaaaaaaand finally, return the lovely path.
}

#pragma endregion

#pragma region Optimizing the Seismic and Pattern Coordinates Paths

vector<Coord>* ToolBox::OptimizeSeismicCoordinates(double*** seismicConditioning, vector<Coord>& simulationPath, ConfigFile* configFile)
{
	return &simulationPath;										// Nothing here to be done (yet!)
}

vector<Coord>* ToolBox::OptimizePatternCoordinates(short*** trainingImage, vector<Coord>& patternsPath, ConfigFile* configFile)
{
	return &patternsPath;										// Nothing here to be done (yet!)
}

#pragma endregion





// Associate the Patterns from the TI to the Seismic Data Events, by comparisons between the SDE and the Synthetic Seismic Data that corresponds to the PTI.
unordered_map<Coord, vector<Coord>*, Coord>* ToolBox::Associate(
	vector<Coord>& seismicDataEventsCoords,					// Coordinates of the Seismic Data Events, inside the Seismic Conditioning matrix,
	double*** seismicConditioning,							// Seismic Conditioning matrix,
	vector<Coord>& patternCoords,							// Coordinates of the Patterns, inside the Training Image AND the Synthetic Seismic Data,
	double*** syntheticSeismicData,							// Synthetic Seismic Data (represented by Optimized Patterns from the Training Image),
	ConfigFile* configFile)									// and the support data.
{


	//return AssociateByPercentageOfDistance(seismicDataEventsCoords, seismicConditioning, patternCoords, syntheticSeismicData, configFile);
	return AssociateByPercentageOfTheTotal(seismicDataEventsCoords, seismicConditioning, patternCoords, syntheticSeismicData, configFile);


}


unordered_map<Coord, vector<Coord>*, Coord>* ToolBox::AssociateByPercentageOfTheTotal(
	vector<Coord>& seismicDataEventsCoords, double*** seismicConditioning,									// the Optimized Seismic Data Events,
	vector<Coord>& patternCoords, double*** syntheticSeismicData,											// the Synthetic Seismic Data (represented by Optimized Patterns from the Training Image),
	ConfigFile* configFile)		// and the support data.
{
	unordered_map<Coord, vector<Coord>*, Coord>* leMap = new unordered_map<Coord, vector<Coord>*, Coord>();	// the Map that will be returned.

	int templateWidth = configFile->GetTemplateWidth();														// the dimensions of the templates
	int templateHeight = configFile->GetTemplateHeight();													// the dimensions of the templates
	int templateDepth = configFile->GetTemplateDepth();														// the dimensions of the templates

	double*** syntheticMatrix = InitializeMatrix_double(templateWidth, templateHeight, templateDepth);		// placeholder for the Synthetic pattern
	double*** seismicMatrix = InitializeMatrix_double(templateWidth, templateHeight, templateDepth);		// placeholder for the Seismic pattern

	Coord sizeOfBigMatrix = Coord(configFile->GetTIWidth(), configFile->GetTIHeight(), configFile->GetTIDepth());	// holds the size of the Seismic Conditioning Matrix
	Coord sizeOfSmallMatrix = Coord(templateWidth, templateHeight, templateDepth);									// holds the size of the Template

	priority_queue<pair<double, int>, vector<pair<double, int>>, LessThanByFirst<int>> pqueue;				// Initialization of the SMALL Priority Queue.
	int maxAmountOfPatternsInTheQueue = configFile->GetMaxAmountOfPatternsInTheQueue();						// Parameter about the size of the Priority Queue.

	for each (Coord seismicDataEvent in seismicDataEventsCoords)											// iterating through all the Seismic Data Events...
	{
		GetSubMatrix(seismicConditioning, sizeOfBigMatrix, seismicDataEvent, seismicMatrix, sizeOfSmallMatrix);		// we get this SDE...


		for (int iPattern = 0; iPattern < patternCoords.size(); iPattern++)									// iterating through all the patterns...
		{
			GetSubMatrix(syntheticSeismicData, sizeOfBigMatrix, patternCoords[iPattern], syntheticMatrix, sizeOfSmallMatrix);		// we get the current pattern

																																	// Calculate the Distance between the Seismic Conditioning Data Event and the Synthetic Seismic Pattern from the Convoluted Training Image.
			double distance = EuclideanDistance(seismicMatrix, syntheticMatrix, templateWidth, templateHeight, templateDepth);

			pqueue.push(make_pair(distance, iPattern));														// put this pair in the queue

			if (pqueue.size() > maxAmountOfPatternsInTheQueue)												// if the PQueue overflows,
				pqueue.pop();																				// remove the one with the largest Distance (most dissimilar)
		}


		(*leMap)[seismicDataEvent] = new vector<Coord>();													// create the new Vector, and then, fill it!
		for (int iPattern = 0; iPattern < maxAmountOfPatternsInTheQueue; iPattern++)						// with every pattern that was in the PQueue.
		{
			(*leMap)[seismicDataEvent]->push_back(patternCoords[pqueue.top().second]);
			pqueue.pop();
		}

		pqueue = priority_queue<pair<double, int>, vector<pair<double, int>>, LessThanByFirst<int>>();		// clear the queue.
	}

	//// Now here we have the full Distance Matrix between all Seismic Conditioning Data Event and the Synthetic Seismic Pattern from the Convoluted Training Image.
	//double normalizedDistanceIJ = 0;																			// Placeholder for the distance value.
	//double thresholdMaxDistance = configFile->GetThresholdMaxDistance();
	//for (int iPattern = 0; iPattern < patternCoords.size(); iPattern++)
	//{
	//	pattern = patternCoords[iPattern];
	//	for (int jSeismic = 0; jSeismic < seismicDataEventsCoords.size(); jSeismic++)
	//	{
	//		//normalizedDistanceIJ = distanceMatrix[iPattern][jSeismic] / maxDistance;
	//		cout << "Distance between Pattern" << iPattern << " and Seismic " << jSeismic << " is " << normalizedDistanceIJ;
	//		if (normalizedDistanceIJ < thresholdMaxDistance)
	//		{
	//			cout << "   SELECTED!!!";
	//			seismicDataEvent = seismicDataEventsCoords[jSeismic];
	//			if ((*leMap)[seismicDataEvent] == nullptr)
	//				(*leMap)[seismicDataEvent] = new vector<Coord>();
	//			(*leMap)[seismicDataEvent]->push_back(pattern);
	//		}
	//		cout << endl;
	//	}
	//}




	DeleteMatrix(syntheticMatrix, templateWidth, templateHeight, templateDepth);
	DeleteMatrix(seismicMatrix, templateWidth, templateHeight, templateDepth);

	return leMap;
}




unordered_map<Coord, vector<Coord>*, Coord>* ToolBox::AssociateByPercentageOfDistance(
	vector<Coord>& seismicDataEventsCoords, double*** seismicConditioning,									// the Optimized Seismic Data Events,
	vector<Coord>& patternCoords, double*** syntheticSeismicData,											// the Synthetic Seismic Data (represented by Optimized Patterns from the Training Image),
	ConfigFile* configFile)		// and the support data.
{
	unordered_map<Coord, vector<Coord>*, Coord>* leMap = new unordered_map<Coord, vector<Coord>*, Coord>();	// the Map that will be returned.

	int templateWidth = configFile->GetTemplateWidth();													// the dimensions of the templates
	int templateHeight = configFile->GetTemplateHeight();													// the dimensions of the templates
	int templateDepth = configFile->GetTemplateDepth();													// the dimensions of the templates


	double*** syntheticMatrix = InitializeMatrix_double(templateWidth, templateHeight, templateDepth);		// placeholder for the Synthetic pattern
	double*** seismicMatrix = InitializeMatrix_double(templateWidth, templateHeight, templateDepth);		// placeholder for the Seismic pattern

																											//double minCompat = configFile->GetMinimumCompatibility();												// Minimum compatibility...   not sure what I'm doing with it.
	double maxDistance = 0;																					// the Maximum Distance, for normalizing the Distance Matrix.

	Coord sizeOfBigMatrix = Coord(configFile->GetTIWidth(), configFile->GetTIHeight(), configFile->GetTIDepth());	// holds the size of the Seismic Conditioning Matrix
	Coord sizeOfSmallMatrix = Coord(templateWidth, templateHeight, templateDepth);									// holds the size of the Template

	double** distanceMatrix = InitializeMatrix_double((int)patternCoords.size(), (int)seismicDataEventsCoords.size());


	for (int iPattern = 0; iPattern < patternCoords.size(); iPattern++)										// iterating through all the patterns...
	{
		GetSubMatrix(syntheticSeismicData, sizeOfBigMatrix, patternCoords[iPattern], syntheticMatrix, sizeOfSmallMatrix);			// we get the current pattern

		for (int jSeismic = 0; jSeismic < seismicDataEventsCoords.size(); jSeismic++)						// iterating through all the Seismic Data Events...
		{
			GetSubMatrix(seismicConditioning, sizeOfBigMatrix, seismicDataEventsCoords[jSeismic], seismicMatrix, sizeOfSmallMatrix);	// we get this SDE...

																																		// Calculate the Distance between the Seismic Conditioning Data Event and the Synthetic Seismic Pattern from the Convoluted Training Image.
			double distance = EuclideanDistance(seismicMatrix, syntheticMatrix, templateWidth, templateHeight, templateDepth);
			//double compatibility = CalculateCompatibility(seismicMatrix, syntheticMatrix, templateWidth, templateHeight, templateDepth, configFile->GetSigmaSquared());

			distanceMatrix[iPattern][jSeismic] = distance;													// store this Distance on the Distance Matrix

			if (distance > maxDistance)																		// and update the maxDistance if needed.
				maxDistance = distance;
		}
	}

	// Now here we have the full Distance Matrix between all Seismic Conditioning Data Event and the Synthetic Seismic Pattern from the Convoluted Training Image.
	Coord seismicDataEvent;																					// The Coords to store data into the Map.
	Coord pattern;																							// The Coords to store data into the Map.
	double normalizedDistanceIJ;																			// Placeholder for the distance value.
	double thresholdMaxDistance = configFile->GetThresholdMaxDistance();

	for (int iPattern = 0; iPattern < patternCoords.size(); iPattern++)
	{
		pattern = patternCoords[iPattern];
		for (int jSeismic = 0; jSeismic < seismicDataEventsCoords.size(); jSeismic++)
		{
			normalizedDistanceIJ = distanceMatrix[iPattern][jSeismic] / maxDistance;

			cout << "Distance between Pattern" << iPattern << " and Seismic " << jSeismic << " is " << normalizedDistanceIJ;

			if (normalizedDistanceIJ < thresholdMaxDistance)
			{
				cout << "   SELECTED!!!";
				seismicDataEvent = seismicDataEventsCoords[jSeismic];

				if ((*leMap)[seismicDataEvent] == nullptr)
					(*leMap)[seismicDataEvent] = new vector<Coord>();

				(*leMap)[seismicDataEvent]->push_back(pattern);
			}

			cout << endl;
		}
	}




	DeleteMatrix(syntheticMatrix, templateWidth, templateHeight, templateDepth);
	DeleteMatrix(seismicMatrix, templateWidth, templateHeight, templateDepth);
	DeleteMatrix(distanceMatrix, (int)patternCoords.size(), (int)seismicDataEventsCoords.size());

	return leMap;
}






void ToolBox::RunSimulationWithSeismicConditioning(
	short*** simulationGrid,												// The Simulation Grid, comes empty, goes full.
	vector<Coord>& simulationPath,											// The Simulation Path, guides the execution of the Simulation.
	unordered_map<Coord, vector<Coord>*, Coord>& association,				// The Association, knows which Patterns correspond to which Data Event.
	double*** seismicConditioning,											// The Seismic Conditioning........................................... I don't know what is this doing here.
	short*** trainingImage,													// The Training Image, to get the Patterns given their Coordinates.
	ConfigFile* configFile)													// And finally, the ConfigFile, for any further configuration.
{
	// Initialize some stuff here....
	int templateWidth = configFile->GetTemplateWidth();
	int templateHeight = configFile->GetTemplateHeight();
	int templateDepth = configFile->GetTemplateDepth();
	Coord templateSize = Coord(templateWidth, templateHeight, templateDepth);
	Coord simGridSize = Coord(configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	Coord tiSize = Coord(configFile->GetTIWidth(), configFile->GetTIHeight(), configFile->GetTIDepth());	// holds the size of the Seismic Conditioning Matrix
	short*** dataEventMatrix = InitializeMatrix_short(templateWidth, templateHeight, templateDepth);
	short*** selectedPattern = InitializeMatrix_short(templateWidth, templateHeight, templateDepth);
	Coord dataEventCoord;

	cout << "We have MEBC_3D!!!!!!!!!!" << endl;


	/// THE FIRST DATA EVENT OF THE SIMULATION PATH IS DONE BY HAND!!!!!!!!!!!!!!
	Coord firstDataEventCoord = simulationPath[0];											// coordenates on the SimGrid of the First Data Event
	vector<Coord> associatedPatterns = *association[firstDataEventCoord];					// retrieving the Patterns associated to the First Data Event

	int randomIndexFromPatternsPath = rand() % associatedPatterns.size();					// from the Patterns of the first Data Event, we select one randomly,
	Coord bestPatternCoord = associatedPatterns[randomIndexFromPatternsPath];				// here it is...
																							//cout << "Selected the pattern " << randomIndexFromPatternsPath << endl;

	GetSubMatrix(trainingImage, tiSize, bestPatternCoord, selectedPattern, templateSize);			// so we 'pull' the Pattern out from the Training Image,
	PasteSubMatrix(simulationGrid, simGridSize, firstDataEventCoord, selectedPattern, templateSize);		// and we 'paste' it into the Simulation Grid.


																											/// THE REST OF THE DATA EVENTS ARE DONE IN HERE!!!!!!!!!!!!!!
	for (int jDataEvent = 1; jDataEvent < simulationPath.size(); jDataEvent++)				// Here comes the main loop of the Simulation, starting from the second DE...
	{
		dataEventCoord = simulationPath[jDataEvent];										// For each Data Event remaining on the Simulation Path...
		GetSubMatrix(simulationGrid, simGridSize, dataEventCoord, dataEventMatrix, templateSize);			// Get the Data Event itself (matrix) from the Simulation Grid,

		associatedPatterns = *association[dataEventCoord];									// Get the Associated Patterns from the TI to this Data Event,

		SelectPattern(trainingImage, associatedPatterns, dataEventMatrix, selectedPattern, configFile);	// Pass the Data Event, Patterns and TI to this Selection method here,

		PasteSubMatrix_MEBC(simulationGrid, simGridSize, dataEventCoord, selectedPattern, templateSize, configFile);	// and finally paste the selected Pattern from the TI into the Simulation Grid.
	}

	DeleteMatrix(dataEventMatrix, templateWidth, templateHeight, templateDepth);			// a little cleanup here...
	DeleteMatrix(selectedPattern, templateWidth, templateHeight, templateDepth);			// a little cleanup there...
}


void ToolBox::RunSimulationWithSeismicConditioning_SeismicImage(
	short*** simulationGrid,												// The Simulation Grid, comes empty, goes full.
	vector<Coord>& simulationPath,											// The Simulation Path, guides the execution of the Simulation.
	unordered_map<Coord, vector<Coord>*, Coord>& association,				// The Association, knows which Patterns correspond to which Data Event.
	double*** seismicConditioning,											// The Seismic Conditioning........................................... I don't know what is this doing here.
	short*** trainingImage,													// The Training Image, to get the Patterns given their Coordinates.
	double*** syntheticSeismicData,											// The Synthetic Seismic Data, to build a "Seismic Realization" parallel to the SimGrid.
	int currentSimulation,
	ConfigFile* configFile)													// And finally, the ConfigFile, for any further configuration.

{
	// Initialize some stuff here....
	int templateWidth = configFile->GetTemplateWidth();
	int templateHeight = configFile->GetTemplateHeight();
	int templateDepth = configFile->GetTemplateDepth();
	Coord templateSize = Coord(templateWidth, templateHeight, templateDepth);
	Coord simGridSize = Coord(configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	Coord tiSize = Coord(configFile->GetTIWidth(), configFile->GetTIHeight(), configFile->GetTIDepth());	// holds the size of the Seismic Conditioning Matrix
	short*** dataEventMatrix = InitializeMatrix_short(templateWidth, templateHeight, templateDepth);
	short*** selectedPattern = InitializeMatrix_short(templateWidth, templateHeight, templateDepth);
	double*** selectedSeismicPattern = InitializeMatrix_double(templateWidth, templateHeight, templateDepth);
	double*** synthSeismicSimGrid = InitializeMatrix_double(configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	Coord dataEventCoord;


	/// THE FIRST DATA EVENT OF THE SIMULATION PATH IS DONE BY HAND!!!!!!!!!!!!!!
	Coord firstDataEventCoord = simulationPath[0];											// coordenates on the SimGrid of the First Data Event
	vector<Coord> associatedPatterns = *association[firstDataEventCoord];					// retrieving the Patterns associated to the First Data Event

	int randomIndexFromPatternsPath = rand() % associatedPatterns.size();					// from the Patterns of the first Data Event, we select one randomly,
	Coord bestPatternCoord = associatedPatterns[randomIndexFromPatternsPath];				// here it is...
	cout << "Selected the pattern " << randomIndexFromPatternsPath << endl;

	GetSubMatrix(trainingImage, tiSize, bestPatternCoord, selectedPattern, templateSize);				// so we 'pull' the Pattern out from the Training Image,
	PasteSubMatrix(simulationGrid, simGridSize, firstDataEventCoord, selectedPattern, templateSize);		// and we 'paste' it into the Simulation Grid.


																											/// THE REST OF THE DATA EVENTS ARE DONE IN HERE!!!!!!!!!!!!!!
	for (int jDataEvent = 1; jDataEvent < simulationPath.size(); jDataEvent++)				// Here comes the main loop of the Simulation, starting from the second DE...
	{
		dataEventCoord = simulationPath[jDataEvent];										// For each Data Event remaining on the Simulation Path...
		GetSubMatrix(simulationGrid, simGridSize, dataEventCoord, dataEventMatrix, templateSize);			// Get the Data Event itself (matrix) from the Simulation Grid,

		associatedPatterns = *association[dataEventCoord];									// Get the Associated Patterns from the TI to this Data Event,

																							//SelectPattern(trainingImage, associatedPatterns, dataEventMatrix, selectedPattern, configFile);	// Pass the Data Event, Patterns and TI to this Selection method here,
		bestPatternCoord = SelectPattern(trainingImage, associatedPatterns, dataEventMatrix, configFile);	// Pass the Data Event, Patterns and TI to this Selection method here,

		GetSubMatrix(trainingImage, tiSize, bestPatternCoord, selectedPattern, templateSize);			// get the selected Pattern from the Best Pat Coord and the TI...
		PasteSubMatrix(simulationGrid, simGridSize, dataEventCoord, selectedPattern, templateSize);		// and finally paste the selected Pattern from the TI into the Simulation Grid.

		GetSubMatrix(syntheticSeismicData, simGridSize, bestPatternCoord, selectedSeismicPattern, templateSize);
		PasteSubMatrix(synthSeismicSimGrid, simGridSize, dataEventCoord, selectedSeismicPattern, templateSize);		// and finally paste the selected Pattern from the TI into the Simulation Grid.
	}

	SaveSimulationToSGEMS(synthSeismicSimGrid, configFile, currentSimulation);
	cout << "	Seismic Simulation saved!" << endl;								// here we have saved it to a new file...

	DeleteMatrix(dataEventMatrix, templateWidth, templateHeight, templateDepth);			// a little cleanup here...
	DeleteMatrix(selectedPattern, templateWidth, templateHeight, templateDepth);			// a little cleanup there...
	DeleteMatrix(selectedSeismicPattern, templateWidth, templateHeight, templateDepth);			// a little cleanup there...
	DeleteMatrix(synthSeismicSimGrid, configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
}


// EXHAUSTIVE SEARCH!!!!!!!!!!!!!!!!!!!!!!!!
void ToolBox::SelectPattern(short*** trainingImage, vector<Coord>& associatedPatterns, short*** dataEventMatrix, short*** selectedPattern, ConfigFile* configFile)
{
	// Initialize some stuff here....
	int templateWidth = configFile->GetTemplateWidth();
	int templateHeight = configFile->GetTemplateHeight();
	int templateDepth = configFile->GetTemplateDepth();
	short*** patternMatrix = InitializeMatrix_short(templateWidth, templateHeight, templateDepth);	// Initialize the placeholder for the current Pattern.
	Coord templateSize = Coord(templateWidth, templateHeight, templateDepth);
	Coord tiSize = Coord(configFile->GetTIWidth(), configFile->GetTIHeight(), configFile->GetTIDepth());	// holds the size of the Seismic Conditioning Matrix

	int maxCandidates = configFile->GetMaxCandidates();
	priority_queue<pair<double, int>, vector<pair<double, int>>, LessThanByFirst<int>> pqueue;		// Initialize the Priority Queue with the largest on Top.

																									// EXHAUSTIVE SEARCH!!!!!!!!!!!!!!!!!!!!!!!!
																									//for each (Coord patternCoord in associatedPatterns)											// Iterating over all given Patterns...
	for (int iPattern = 0; iPattern < associatedPatterns.size(); iPattern++)
	{
		GetSubMatrix(trainingImage, tiSize, associatedPatterns[iPattern], patternMatrix, templateSize);		// Get the current Pattern,

																											// Now, compare the given Data Event to the current Pattern
		double distance = EuclideanDistance_WithoutZeroes(dataEventMatrix, patternMatrix, templateWidth, templateHeight, templateDepth);

		pqueue.push(make_pair(distance, iPattern));													// put this <Distance, Coordinate> pair on the PQueue,

		if (pqueue.size() > maxCandidates)															// if the PQueue overflows,
			pqueue.pop();																			// remove the one with the largest Distance (most dissimilar)
	}


	int selectedIndex = rand() % maxCandidates;														// randomly select one of the most similar
	for (int i = 0; i < selectedIndex; i++)															// advance the PQueue to it...
		pqueue.pop();

	Coord bestPatternCoord = associatedPatterns[pqueue.top().second];								// and mark it as the Best Pattern.
																									//cout << "Selected the pattern " << pqueue.top().second << endl;


																									// Now here we have the bestPatternCoord, it's just getting the Pattern from the TI into the selectedPattern
	GetSubMatrix(trainingImage, tiSize, bestPatternCoord, selectedPattern, templateSize);

	DeleteMatrix(patternMatrix, templateWidth, templateHeight, templateDepth);
}


Coord ToolBox::SelectPattern(short*** trainingImage, vector<Coord>& associatedPatterns, short*** dataEventMatrix, ConfigFile* configFile)
{
	// Initialize some stuff here....
	int templateWidth = configFile->GetTemplateWidth();
	int templateHeight = configFile->GetTemplateHeight();
	int templateDepth = configFile->GetTemplateDepth();
	short*** patternMatrix = InitializeMatrix_short(templateWidth, templateHeight, templateDepth);	// Initialize the placeholder for the current Pattern.
	Coord templateSize = Coord(templateWidth, templateHeight, templateDepth);
	Coord tiSize = Coord(configFile->GetTIWidth(), configFile->GetTIHeight(), configFile->GetTIDepth());	// holds the size of the Seismic Conditioning Matrix

	int maxCandidates = configFile->GetMaxCandidates();
	priority_queue<pair<double, int>, vector<pair<double, int>>, LessThanByFirst<int>> pqueue;		// Initialize the Priority Queue with the largest on Top.

																									// EXHAUSTIVE SEARCH!!!!!!!!!!!!!!!!!!!!!!!!
																									//for each (Coord patternCoord in associatedPatterns)											// Iterating over all given Patterns...
	for (int iPattern = 0; iPattern < associatedPatterns.size(); iPattern++)
	{
		GetSubMatrix(trainingImage, tiSize, associatedPatterns[iPattern], patternMatrix, templateSize);		// Get the current Pattern,

																											// Now, compare the given Data Event to the current Pattern
		double distance = EuclideanDistance_WithoutZeroes(dataEventMatrix, patternMatrix, templateWidth, templateHeight, templateDepth);

		pqueue.push(make_pair(distance, iPattern));													// put this <Distance, Coordinate> pair on the PQueue,

		if (pqueue.size() > maxCandidates)															// if the PQueue overflows,
			pqueue.pop();																			// remove the one with the largest Distance (most dissimilar)
	}


	int selectedIndex = rand() % maxCandidates;														// randomly select one of the most similar
	for (int i = 0; i < selectedIndex; i++)															// advance the PQueue to it...
		pqueue.pop();

	Coord bestPatternCoord = associatedPatterns[pqueue.top().second];								// and mark it as the Best Pattern.
																									//cout << "Selected the pattern " << pqueue.top().second << endl;


																									// Now here we have the bestPatternCoord, it's just getting the Pattern from the TI into the selectedPattern
																									//GetSubMatrix(trainingImage, bestPatternCoord, selectedPattern, configFile);

	DeleteMatrix(patternMatrix, templateWidth, templateHeight, templateDepth);

	return bestPatternCoord;
}


template <class T>
void ToolBox::GetSubMatrix(T*** bigMatrix, Coord& sizeOfBigMatrix, Coord& startInBigMatrix, T*** smallMatrix, Coord& sizeOfSmallMatrix)
{
	int x = startInBigMatrix.x;
	int y = startInBigMatrix.y;
	int z = startInBigMatrix.z;

	int matrixWidth = min(sizeOfSmallMatrix.x, sizeOfBigMatrix.x - x);
	int matrixHeight = min(sizeOfSmallMatrix.y, sizeOfBigMatrix.y - y);
	int matrixDepth = min(sizeOfSmallMatrix.z, sizeOfBigMatrix.z - z);

	for (int d = 0; d < matrixDepth; d++)
		for (int h = 0; h < matrixHeight; h++)
			for (int w = 0; w < matrixWidth; w++)
				smallMatrix[w][h][d] = bigMatrix[x + w][y + h][z + d];

}

template <class T>
void ToolBox::PasteSubMatrix(T*** bigMatrix, Coord& sizeOfBigMatrix, Coord& startInBigMatrix, T*** smallMatrix, Coord& sizeOfSmallMatrix)
{
	int x = startInBigMatrix.x;
	int y = startInBigMatrix.y;
	int z = startInBigMatrix.z;

	int matrixWidth = min(sizeOfSmallMatrix.x, sizeOfBigMatrix.x - x);
	int matrixHeight = min(sizeOfSmallMatrix.y, sizeOfBigMatrix.y - y);
	int matrixDepth = min(sizeOfSmallMatrix.z, sizeOfBigMatrix.z - z);

	for (int d = 0; d < matrixDepth; d++)
		for (int h = 0; h < matrixHeight; h++)
			for (int w = 0; w < matrixWidth; w++)
				bigMatrix[x + w][y + h][z + d] = smallMatrix[w][h][d];

}



//template <class T>
//void ToolBox::PasteSubMatrix_MEBC(T*** bigMatrix, Coord startInBigMatrix, T*** smallMatrix, ConfigFile* configFile)
void ToolBox::PasteSubMatrix_MEBC(short*** bigMatrix, Coord& sizeOfBigMatrix, Coord& startInBigMatrix, short*** smallMatrix, Coord& sizeOfSmallMatrix, ConfigFile* configFile)
{
	short*** AMatrix;											// Declaring the A Matrix.
	short*** BMatrix;											// Declaring the B Matrix.
	short*** answer;											// Declaring the Answer Matrix.
	short*** tempMatrix;										// Declaring the Temporary Matrix.

	int sizeX_template = sizeOfSmallMatrix.x;					// Initializing the size of the Template.
	int sizeY_template = sizeOfSmallMatrix.y;					// Initializing the size of the Template.
	int sizeZ_template = sizeOfSmallMatrix.z;					// Initializing the size of the Template.

	int sizeX_overlap = configFile->GetOverlapWidth();			// Initializing the size of the Overlap.
	int sizeY_overlap = configFile->GetOverlapHeight();			// Initializing the size of the Overlap.
	int sizeZ_overlap = configFile->GetOverlapDepth();			// Initializing the size of the Overlap.

	Coord zeroZeroZero = Coord(0, 0, 0);							// The starting Coordinate in the SmallMatrix.



																	// STEP 1. TOP:		(TxTxO)
	Coord topMatrixDimensions = Coord(sizeX_template, sizeY_template, sizeZ_overlap);					// Initializing the Top Matrix's Dimensions.

	AMatrix = InitializeMatrix_short(topMatrixDimensions);
	GetSubMatrix(bigMatrix, sizeOfBigMatrix, startInBigMatrix, AMatrix, topMatrixDimensions);			// Put into AMatrix the Original's slice.

	BMatrix = InitializeMatrix_short(topMatrixDimensions);
	GetSubMatrix(smallMatrix, sizeOfSmallMatrix, zeroZeroZero, BMatrix, topMatrixDimensions);			// Put into BMatrix the Pattern's slice.

	answer = InitializeMatrix_short(topMatrixDimensions);
	MEBC_3D(AMatrix, BMatrix, sizeX_template, sizeY_template, sizeZ_overlap, answer);					// Calculate the MEBC_3D and put it into Answer Matrix.

	PasteSubMatrix(bigMatrix, sizeOfBigMatrix, startInBigMatrix, answer, topMatrixDimensions);			// Finally, paste the AnswerMatrix into the Original Matrix.






																										// STEP 2. SIDE:	(OxTxT)
	Coord sideMatrixDimensions = Coord(sizeX_overlap, sizeY_template, sizeZ_template);					// Initializing the Side Matrix's Dimensions.
	tempMatrix = InitializeMatrix_short(sideMatrixDimensions);

	GetSubMatrix(bigMatrix, sizeOfBigMatrix, startInBigMatrix, tempMatrix, sideMatrixDimensions);		// Put into TempMatrix the Original's slice.
	FLIP_SIDEMATRIX_INTO_TOPMATRIX(tempMatrix, sideMatrixDimensions, AMatrix, topMatrixDimensions);		// Flip the TempMatrix into AMatrix.

	GetSubMatrix(smallMatrix, sizeOfSmallMatrix, zeroZeroZero, tempMatrix, sideMatrixDimensions);		// Put into TempMatrix the Pattern's slice.
	FLIP_SIDEMATRIX_INTO_TOPMATRIX(tempMatrix, sideMatrixDimensions, BMatrix, topMatrixDimensions);		// Flip the TempMatrix into BMatrix.

	MEBC_3D(AMatrix, BMatrix, sizeX_template, sizeY_template, sizeZ_overlap, answer);					// Calculate the MEBC_3D and put it into Answer Matrix.
	FLIP_TOPMATRIX_INTO_SIDEMATRIX(answer, topMatrixDimensions, tempMatrix, sideMatrixDimensions);		// Flip the AnswerMatrix into TempMatrix.

	PasteSubMatrix(bigMatrix, sizeOfBigMatrix, startInBigMatrix, answer, topMatrixDimensions);			// Finally, paste the TempMatrix into the Original Matrix.
	DeleteMatrix(tempMatrix, sideMatrixDimensions);														// aaaaaand some housekeeping for the TempMatrix.





																										// STEP 3. FRONT:	(TxOxT)
	Coord frontMatrixDimensions = Coord(sizeX_template, sizeY_overlap, sizeZ_template);					// Initializing the front Matrix's Dimensions.
	tempMatrix = InitializeMatrix_short(frontMatrixDimensions);

	GetSubMatrix(bigMatrix, sizeOfBigMatrix, startInBigMatrix, tempMatrix, frontMatrixDimensions);		// Put into TempMatrix the Original's slice.
	FLIP_FRONTMATRIX_INTO_TOPMATRIX(tempMatrix, frontMatrixDimensions, AMatrix, topMatrixDimensions);	// Flip the TempMatrix into AMatrix.

	GetSubMatrix(smallMatrix, sizeOfSmallMatrix, zeroZeroZero, tempMatrix, frontMatrixDimensions);		// Put into TempMatrix the Pattern's slice.
	FLIP_FRONTMATRIX_INTO_TOPMATRIX(tempMatrix, frontMatrixDimensions, BMatrix, topMatrixDimensions);	// Flip the TempMatrix into BMatrix.

	MEBC_3D(AMatrix, BMatrix, sizeX_template, sizeY_template, sizeZ_overlap, answer);					// Calculate the MEBC_3D and put it into Answer Matrix.
	FLIP_TOPMATRIX_INTO_FRONTMATRIX(answer, topMatrixDimensions, tempMatrix, frontMatrixDimensions);	// Flip the AnswerMatrix into TempMatrix.

	PasteSubMatrix(bigMatrix, sizeOfBigMatrix, startInBigMatrix, answer, topMatrixDimensions);			// Finally, paste the TempMatrix into the Original Matrix.
	DeleteMatrix(tempMatrix, frontMatrixDimensions);													// aaaaaand some housekeeping for the TempMatrix.



																										// STEP 4. INSIDE:										FIX IT!!!!!!!!!!!!!!!!!

																										////////////////////////////////////////////////////////////////
	int xPlusO = startInBigMatrix.x + sizeX_overlap;
	int yPlusO = startInBigMatrix.y + sizeY_overlap;
	int zPlusO = startInBigMatrix.z + sizeZ_overlap;
	int matrixWidth = min(sizeOfSmallMatrix.x, sizeOfBigMatrix.x - xPlusO);
	int matrixHeight = min(sizeOfSmallMatrix.y, sizeOfBigMatrix.y - yPlusO);
	int matrixDepth = min(sizeOfSmallMatrix.z, sizeOfBigMatrix.z - zPlusO);
	for (int d = sizeZ_overlap; d < matrixDepth; d++)
		for (int h = sizeY_overlap; h < matrixHeight; h++)
			for (int w = sizeX_overlap; w < matrixWidth; w++)
				bigMatrix[xPlusO + w][yPlusO + h][zPlusO + d] = smallMatrix[w][h][d];
	////////////////////////////////////////////////////////////////




	DeleteMatrix(AMatrix, sizeX_template, sizeY_template, sizeZ_overlap);								// Housekeeping time!!!
	DeleteMatrix(BMatrix, sizeX_template, sizeY_template, sizeZ_overlap);								// Housekeeping time!!!
	DeleteMatrix(answer, sizeX_template, sizeY_template, sizeZ_overlap);								// Housekeeping time!!!
}


void ToolBox::FLIP_SIDEMATRIX_INTO_TOPMATRIX(short*** sideMatrix, Coord& sideMatrixDimensions, short*** topMatrix, Coord& topMatrixDimensions)
{
	for (int x = 0; x < topMatrixDimensions.x; x++)
		for (int y = 0; y < topMatrixDimensions.y; y++)
			for (int z = 0; z < topMatrixDimensions.z; z++)
				topMatrix[x][y][z] = sideMatrix[sideMatrixDimensions.z - z - 1][y][x];
}
void ToolBox::FLIP_TOPMATRIX_INTO_SIDEMATRIX(short*** topMatrix, Coord& topMatrixDimensions, short*** sideMatrix, Coord& sideMatrixDimensions)
{
	for (int x = 0; x < sideMatrixDimensions.x; x++)
		for (int y = 0; y < sideMatrixDimensions.y; y++)
			for (int z = 0; z < sideMatrixDimensions.z; z++)
				sideMatrix[sideMatrixDimensions.z - z - 1][y][x] = topMatrix[x][y][z];
}
void ToolBox::FLIP_FRONTMATRIX_INTO_TOPMATRIX(short*** frontMatrix, Coord& frontMatrixDimensions, short*** topMatrix, Coord& topMatrixDimensions)
{
	for (int x = 0; x < topMatrixDimensions.x; x++)
		for (int y = 0; y < topMatrixDimensions.y; y++)
			for (int z = 0; z < topMatrixDimensions.z; z++)
				topMatrix[x][y][z] = frontMatrix[x][frontMatrixDimensions.z - z - 1][y];
}
void ToolBox::FLIP_TOPMATRIX_INTO_FRONTMATRIX(short*** topMatrix, Coord& topMatrixDimensions, short*** frontMatrix, Coord& frontMatrixDimensions)
{
	for (int x = 0; x < frontMatrixDimensions.x; x++)
		for (int y = 0; y < frontMatrixDimensions.y; y++)
			for (int z = 0; z < frontMatrixDimensions.z; z++)
				frontMatrix[x][frontMatrixDimensions.z - z - 1][y] = topMatrix[x][y][z];
}


template <class T>
void ToolBox::PasteSubMatrix_AveragedOverlap(T*** bigMatrix, Coord startInBigMatrix, T*** smallMatrix, ConfigFile* configFile)
{
	int matrixWidth = configFile->GetTemplateWidth();
	int matrixHeight = configFile->GetTemplateHeight();
	int matrixDepth = configFile->GetTemplateDepth();

	int x = startInBigMatrix.x;
	int y = startInBigMatrix.y;
	int z = startInBigMatrix.z;

	for (int d = 0; d < matrixDepth; d++)
		for (int h = 0; h < matrixHeight; h++)
			for (int w = 0; w < matrixWidth; w++)
				bigMatrix[x + w][y + h][z + d] = smallMatrix[w][h][d];

}




//template <class T>
//T** ToolBox::MEBC_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, T** answer)// , ConfigFile* configFile)
//{
//	// STEP 1. Define an error surface, having a rectangular size TEMPLATESIZE by OVERLAP.
//	double** surfaceError = InitializeMatrix_double(sizeX_template, sizeY_overlap);		// Create the surface error
//	for (int i = 0; i < sizeX_template; i++)											// and initialize it
//		for (int j = 0; j < sizeY_overlap; j++)
//			surfaceError[i][j] = (AMatrix[i][j] - BMatrix[i][j])*(AMatrix[i][j] - BMatrix[i][j]); // excellent initialization.
//
//
//	// STEP 2. Compute the cumulative minimum error along the cutting direction
//	double** cumulMinError = InitializeMatrix_double(sizeX_template, sizeY_overlap);	// Create the Cumulative Minimum Error table.
//	double leMinimum = 0;
//
//	for (int i = 0; i < sizeX_template; i++)
//	{
//		for (int j = 0; j < sizeY_overlap; j++)
//			if (i == 0)
//				cumulMinError[0][j] = surfaceError[i][j];
//			else
//			{
//				leMinimum = cumulMinError[i - 1][j];								// initialize as the value in the previous row, same column
//				if (j > 0)
//					leMinimum = min(leMinimum, cumulMinError[i - 1][j - 1]);		// if there's a anterior column, consider it.
//				if (j < sizeY_overlap - 1)
//					leMinimum = min(leMinimum, cumulMinError[i - 1][j + 1]);		// if there's a posterior column, consider it.
//
//				cumulMinError[i][j] = surfaceError[i][j] + leMinimum;			// aaaaand finally, the Cumulative Minimum Error is the current value + the minimum.
//			}
//	}
//
//
//	// STEP 3. Identify the coordinate J corresponding to the entry with smallest value on the last row of E. This location	corresponds to the arrival point of a path of minimum cost through the error surface.
//	int currentI = sizeX_template - 1;										// the analysis starts in the last row of the Cumulative Minimum Error table.
//
//	int * boundary = new int[sizeX_template];								// defining the *boundary*
//	int smallestJIndex = 0;													// and this guy will hold the column index of the smallest value in the current row.
//
//	for (int j = 1; j < sizeY_overlap; j++)									// we iterate over the columns of the current row...
//		if (cumulMinError[currentI][j] < cumulMinError[currentI][smallestJIndex])	// and if there's a column holding a smaller value than the smallest known...
//			smallestJIndex = j;														// Improvise. Adapt. Overcome.  Or simply, memorize the new smallest.
//
//	boundary[currentI] = smallestJIndex;									// and don't forget to put the smallest guy in the last position of the boundary!!!
//
//
//
//	// STEP 4. Trace back the minimum values for each row i, going backward (with i=p-1...1) and each time identify the cutting path.
//	int previousJIndex;														// now, this 'temp' guy will facilitate the readability of the code.
//
//	for (currentI--; currentI >= 0; currentI--)								// iterating over the rest of the rows, in reverse order...
//	{
//		previousJIndex = boundary[currentI + 1];							// get the column index of the   previous row smallest value
//		smallestJIndex = boundary[currentI + 1];							// and consider it as a possible current  row smallest value
//
//																			// if the smallest of the previous row had a column to the 'left' of it, and if the value in this 'left' column is smaller than the smallest known...
//		if (previousJIndex > 0 && cumulMinError[currentI][previousJIndex - 1] < cumulMinError[currentI][smallestJIndex])
//			smallestJIndex = previousJIndex - 1;
//
//		// if the smallest of the previous row had a column to the 'right' of it, and if the value in this 'right' column is smaller than the smallest known...
//		if (previousJIndex < sizeY_overlap - 1 && cumulMinError[currentI][previousJIndex + 1] < cumulMinError[currentI][smallestJIndex])
//			smallestJIndex = previousJIndex + 1;
//
//		boundary[currentI] = smallestJIndex;								// Include the column index of the smallest value of this row in the boundary.
//	}
//
//
//	// STEP 5. Merge the two Surfaces (Matrices) using the calculated Boundary.
//	for (int i = 0; i < sizeX_template; i++)
//	{
//		for (int j = 0; j < boundary[i]; j++)
//			answer[i][j] = AMatrix[i][j];
//		for (int j = boundary[i]; j < sizeY_overlap; j++)
//			answer[i][j] = BMatrix[i][j];
//	}
//
//
//
//	DeleteMatrix(surfaceError, sizeX_template, sizeY_overlap);
//	DeleteMatrix(cumulMinError, sizeX_template, sizeY_overlap);
//	delete[] boundary;
//
//
//	return answer;
//}

//template <class T>
//void ToolBox::PasteSubMatrix_MEBC_Mahmud(T*** bigMatrix, Coord startInBigMatrix, T*** smallMatrix, ConfigFile* configFile)
//{
//	int matrixWidth = configFile->GetTemplateWidth();
//	int matrixHeight = configFile->GetTemplateHeight();
//	int matrixDepth = configFile->GetTemplateDepth();
//
//	int x = startInBigMatrix.x;
//	int y = startInBigMatrix.y;
//	int z = startInBigMatrix.z;
//
//	int olW = configFile->GetOverlapWidth();
//	int olH = configFile->GetOverlapHeight();
//	int olD = configFile->GetOverlapDepth();
//
//
//	// First, paste the non-conflicted area (green area)
//	for (int d = olD; d < matrixDepth; d++)
//		for (int h = olH; h < matrixHeight; h++)
//			for (int w = olW; w < matrixWidth; w++)
//				bigMatrix[x + w][y + h][z + d] = smallMatrix[w][h][d];
//
//
//	// Second, get each of the overlap regions and run the MEBC on them
//
//
//
//	// and Third, paste the 'merged' regions into the bigMatrix
//
//
//
//
//	for (int d = 0; d < matrixDepth; d++)
//		for (int h = 0; h < matrixHeight; h++)
//			for (int w = 0; w < matrixWidth; w++)
//				bigMatrix[x + w][y + h][z + d] = smallMatrix[w][h][d];
//
//}



//template <class T>
//void ToolBox::PasteSubMatrix_MEBC(T*** bigMatrix, Coord startInBigMatrix, T*** smallMatrix, ConfigFile* configFile)
//{
//	int matrixWidth = configFile->GetTemplateWidth();
//	int matrixHeight = configFile->GetTemplateHeight();
//	int matrixDepth = configFile->GetTemplateDepth();
//
//	int x = startInBigMatrix.x;
//	int y = startInBigMatrix.y;
//	int z = startInBigMatrix.z;
//
//	for (int d = 0; d < matrixDepth; d++)
//		for (int h = 0; h < matrixHeight; h++)
//			for (int w = 0; w < matrixWidth; w++)
//				bigMatrix[x + w][y + h][z + d] = smallMatrix[w][h][d];
//
//}


void ToolBox::SaveSimulationToSGEMS(short*** simulationGrid, ConfigFile* configFile, int currentSimulation)
{
	string arqName = "Realizations3D";
	Create_Directory("Realizations3D");

	arqName = arqName + "\\Simulation_" + to_string(currentSimulation) + "_" //+ configFile->GetRealizationFileName()
		+ "R" + to_string(configFile->GetRealizationWidth()) + "x" + to_string(configFile->GetRealizationHeight()) + "x" + to_string(configFile->GetRealizationDepth()) + "_"
		+ "T" + to_string(configFile->GetTemplateWidth()) + "x" + to_string(configFile->GetTemplateHeight()) + "x" + to_string(configFile->GetTemplateDepth()) + "_"
		+ "O" + to_string(configFile->GetOverlapWidth()) + "x" + to_string(configFile->GetOverlapHeight()) + "x" + to_string(configFile->GetOverlapDepth()) //   +"_"
		+ ".sgems";


	//FILE* realizationFile = fopen(arqName.c_str(), "w");
	ofstream myfile;
	myfile.open(arqName);

	//fprintf(realizationFile, "%d %d %d\n", configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	myfile << configFile->GetRealizationWidth() << " " << configFile->GetRealizationHeight() << " " << configFile->GetRealizationDepth() << endl;
	//fprintf(realizationFile, "1\nfacies\n");
	myfile << "1" << endl << "nfacies" << endl;

	for (int d = 0; d < configFile->GetRealizationDepth(); d++)
		for (int h = 0; h < configFile->GetRealizationHeight(); h++)
			for (int w = 0; w < configFile->GetRealizationWidth(); w++)
			{
				//fprintf(realizationFile, "%d\n", simulationGrid[w][h][d]);
				myfile << simulationGrid[w][configFile->GetRealizationHeight() - 1 - h][d] << endl; // Notice that in sgems we need to invert the Y axis (for NO reason)
			}

	//fclose(realizationFile);
	myfile.close();

}


void ToolBox::SaveSimulationToSGEMS(double*** simulationGrid, ConfigFile* configFile, int currentSimulation)
{
	string arqName = "Realizations3D\\";

	arqName = arqName + "SimulationCont_" + to_string(currentSimulation) + "_" //+ configFile->GetRealizationFileName()
		+ "R" + to_string(configFile->GetRealizationWidth()) + "x" + to_string(configFile->GetRealizationHeight()) + "x" + to_string(configFile->GetRealizationDepth()) + "_"
		+ "T" + to_string(configFile->GetTemplateWidth()) + "x" + to_string(configFile->GetTemplateHeight()) + "x" + to_string(configFile->GetTemplateDepth()) + "_"
		+ "O" + to_string(configFile->GetOverlapWidth()) + "x" + to_string(configFile->GetOverlapHeight()) + "x" + to_string(configFile->GetOverlapDepth()) //   +"_"
		+ ".sgems";


	//FILE* realizationFile = fopen(arqName.c_str(), "w");
	ofstream myfile;
	myfile.open(arqName);

	//fprintf(realizationFile, "%d %d %d\n", configFile->GetRealizationWidth(), configFile->GetRealizationHeight(), configFile->GetRealizationDepth());
	myfile << configFile->GetRealizationWidth() << " " << configFile->GetRealizationHeight() << " " << configFile->GetRealizationDepth() << endl;
	//fprintf(realizationFile, "1\nfacies\n");
	myfile << "1" << endl << "nfacies" << endl;

	for (int d = 0; d < configFile->GetRealizationDepth(); d++)
		for (int h = 0; h < configFile->GetRealizationHeight(); h++)
			for (int w = 0; w < configFile->GetRealizationWidth(); w++)
			{
				//fprintf(realizationFile, "%d\n", simulationGrid[w][h][d]);
				myfile << simulationGrid[w][h][d] << endl;
			}

	//fclose(realizationFile);
	myfile.close();

}


#pragma region Compatibility and Distance functions

double ToolBox::EuclideanDistance(double*** matrixA, double*** matrixB, int width, int height, int depth)
{
	double ans = 0;

	for (int d = 0; d < depth; d++)
	{
		for (int h = 0; h < height; h++)
		{
			for (int w = 0; w < width; w++)
			{
				double a = matrixA[w][h][d];
				double b = matrixB[w][h][d];
				double c = a - b;
				ans += pow(a - b, 2);
			}
		}
	}

	ans = sqrt(ans);

	return ans;
}

double ToolBox::EuclideanDistance(short*** matrixA, short*** matrixB, int width, int height, int depth)
{
	double ans = 0;

	for (int d = 0; d < depth; d++)
	{
		for (int h = 0; h < height; h++)
		{
			for (int w = 0; w < width; w++)
			{
				short a = matrixA[w][h][d];
				short b = matrixB[w][h][d];
				short c = a - b;
				ans += pow(a - b, 2);
			}
		}
	}

	ans = sqrt(ans);

	return ans;
}

double ToolBox::EuclideanDistance_WithoutZeroes(short*** matrixA, short*** matrixB, int width, int height, int depth)
{
	double ans = 0;

	for (int d = 0; d < depth; d++)
	{
		for (int h = 0; h < height; h++)
		{
			for (int w = 0; w < width; w++)
			{
				short a = matrixA[w][h][d];
				if (a == 0) continue;

				short b = matrixB[w][h][d];
				if (b == 0) continue;

				short c = a - b;
				ans += pow(a - b, 2);
			}
		}
	}

	ans = sqrt(ans);

	return ans;
}


double ToolBox::CalculateCompatibility(double*** matrixA, double*** matrixB, int width, int height, int depth, double sigmaSquared)
{
	double ans = 0;

	for (int d = 0; d < depth; d++)
	{
		for (int h = 0; h < height; h++)
		{
			for (int w = 0; w < width; w++)
			{
				double a = matrixA[w][h][d];
				double b = matrixB[w][h][d];
				double c = a - b;
				ans += pow(a - b, 2);
			}
		}
	}

	ans = -ans / (2 * sigmaSquared);

	ans = exp(ans);

	return ans;
}

#pragma endregion

#pragma region Initializing, Deleting and Resetting Matrices

short*** ToolBox::InitializeMatrix_short(int width, int height, int depth)
{
	short*** matrix = new short**[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new short*[height];
		for (int j = 0; j < height; j++)
		{
			matrix[i][j] = new short[depth];
			for (int k = 0; k < depth; k++)
				matrix[i][j][k] = 0;
		}
	}

	return matrix;
}

short*** ToolBox::InitializeMatrix_short(Coord& dimensions)
{
	short*** matrix = new short**[dimensions.x];
	for (int i = 0; i < dimensions.x; i++)
	{
		matrix[i] = new short*[dimensions.y];
		for (int j = 0; j < dimensions.y; j++)
		{
			matrix[i][j] = new short[dimensions.z];
			for (int k = 0; k < dimensions.z; k++)
				matrix[i][j][k] = 0;
		}
	}

	return matrix;
}


short** ToolBox::InitializeMatrix_short(int width, int height)
{
	short** matrix = new short*[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new short[height];
		for (int j = 0; j < height; j++)
			matrix[i][j] = 0;
	}

	return matrix;
}

int ** ToolBox::InitializeMatrix_int(int width, int height)
{
	int** matrix = new int*[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new int[height];
		for (int j = 0; j < height; j++)
			matrix[i][j] = 0;
	}

	return matrix;
}


short*** ToolBox::InitializeMatrix_short_Random(int width, int height, int depth, int minRandom, int maxRandom)
{
	short*** matrix = new short**[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new short*[height];
		for (int j = 0; j < height; j++)
		{
			matrix[i][j] = new short[depth];
			for (int k = 0; k < depth; k++)
				matrix[i][j][k] = rand() % (maxRandom - minRandom + 1) + minRandom;
		}
	}

	return matrix;
}

short** ToolBox::InitializeMatrix_short_Random(int width, int height, int minRandom, int maxRandom)
{
	short** matrix = new short*[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new short[height];
		for (int j = 0; j < height; j++)
			matrix[i][j] = rand() % (maxRandom - minRandom + 1) + minRandom;
	}

	return matrix;
}


double*** ToolBox::InitializeMatrix_double(int width, int height, int depth)
{
	double*** matrix = new double**[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new double*[height];
		for (int j = 0; j < height; j++)
		{
			matrix[i][j] = new double[depth];
			for (int k = 0; k < depth; k++)
				matrix[i][j][k] = 0;
		}
	}

	return matrix;
}

double** ToolBox::InitializeMatrix_double(int width, int height)
{
	double** matrix = new double*[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new double[height];
		for (int j = 0; j < height; j++)
			matrix[i][j] = 0;
	}

	return matrix;
}


double*** ToolBox::InitializeMatrix_double_Random(int width, int height, int depth, int minRandom, int maxRandom)
{
	double*** matrix = new double**[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new double*[height];
		for (int j = 0; j < height; j++)
		{
			matrix[i][j] = new double[depth];
			for (int k = 0; k < depth; k++)
				matrix[i][j][k] = rand() % (maxRandom - minRandom + 1) + minRandom;
		}
	}

	return matrix;
}

double** ToolBox::InitializeMatrix_double_Random(int width, int height, int minRandom, int maxRandom)
{
	double** matrix = new double*[width];
	for (int i = 0; i < width; i++)
	{
		matrix[i] = new double[height];
		for (int j = 0; j < height; j++)
			matrix[i][j] = rand() % (maxRandom - minRandom + 1) + minRandom;
	}

	return matrix;
}




//template <class T>
//void ToolBox::DeleteMatrix(T*** matrix, int width, int height, int depth)
//{
//	for (int i = 0; i < width; i++)
//	{
//		for (int j = 0; j < height; j++)
//		{
//			//for (int k = 0; k < depth; k++)
//			//{
//			//	delete matrix[i][j][k];
//			//}
//			delete[] matrix[i][j];
//
//		}
//		delete[] matrix[i];
//	}
//
//	delete[] matrix;
//}
void ToolBox::DeleteMatrix(short*** matrix, Coord& dimensions)
{
	for (int i = 0; i < dimensions.x; i++)
	{
		for (int j = 0; j < dimensions.y; j++)
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


void ToolBox::DeleteMatrix(short*** matrix, int width, int height, int depth)
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

void ToolBox::DeleteMatrix(double*** matrix, int width, int height, int depth)
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

void ToolBox::DeleteMatrix(double** matrix, int width, int height)
{
	for (int i = 0; i < width; i++)
	{
		//for (int j = 0; j < height; j++)
		//{
		//	//for (int k = 0; k < depth; k++)
		//	//{
		//	//	delete matrix[i][j][k];
		//	//}
		//	delete[] matrix[i][j];

		//}
		delete[] matrix[i];
	}

	delete[] matrix;
}

void ToolBox::DeleteMatrix(int** matrix, int width, int height)
{
	for (int i = 0; i < width; i++)
	{
		//for (int j = 0; j < height; j++)
		//{
		//	//for (int k = 0; k < depth; k++)
		//	//{
		//	//	delete matrix[i][j][k];
		//	//}
		//	delete[] matrix[i][j];

		//}
		delete[] matrix[i];
	}

	delete[] matrix;
}


template <class T>
void ToolBox::ResetMatrix(T*** matrix, int width, int height, int depth)
{
	for (int d = 0; d < depth; d++)
		for (int h = 0; h < height; h++)
			for (int w = 0; w < width; w++)
				matrix[w][h][d] = 0;

}
void ToolBox::ResetMatrix(short*** matrix, int width, int height, int depth)
{
	for (int d = 0; d < depth; d++)
		for (int h = 0; h < height; h++)
			for (int w = 0; w < width; w++)
				matrix[w][h][d] = 0;

}
//
//void ToolBox::ResetMatrix(double*** matrix, int width, int height, int depth)
//{
//	for (int d = 0; d < depth; d++)
//		for (int h = 0; h < height; h++)
//			for (int w = 0; w < width; w++)
//				matrix[w][h][d] = 0;
//
//}

#pragma endregion

#pragma region PAIRING AND DEPAIRING, ACCORDING TO GEORG CANTOR

int ToolBox::Pair_CantorPairing(int x, int y)
{
	int w = x + y;
	int pairr = w * (w + 1) / 2 + y;
	return pairr;
}
pair<int, int>* ToolBox::Depair_CantorPairing(int z)
{
	double aux1 = 8;
	aux1 = aux1 * z;
	aux1 = aux1 + 1;
	aux1 = sqrt(aux1);
	aux1 = aux1 - 1;
	aux1 = aux1 / 2;

	int w = (int)aux1;

	int t = w * (w + 1) / 2;

	int y = z - t;
	int x = w - y;


	pair<int, int>* couple = new pair<int, int>(x, y);
	return couple;
}
void ToolBox::TestPairingDepairing(int maxLimit)
{
	int pairedValues = 0, x = 0, y = 0;
	pair<int, int>* depairedValues;

	for (int i = 0; i < maxLimit; i++)
	{
		for (int j = 0; j < maxLimit; j++)
		{
			pairedValues = Pair_CantorPairing(i, j);
			depairedValues = Depair_CantorPairing(pairedValues);

			if (i != depairedValues->first || j != depairedValues->second)
				cout << "WE HAVE AN ERROR WITH THE PAIRING OF " << i << "," << j << endl;

		}
	}
}

bool ToolBox::Create_Directory(char* folderName)
{
	const char *fn = folderName;								 // Converting to a type readable for CreateDirectory function
	size_t size = strlen(fn) + 1;								 // ...
	wchar_t *folder = new wchar_t[size]();						 // ...
	size_t outSize;												 // ...
	mbstowcs_s(&outSize, folder, size, fn, size - 1);			 // finishing converting

	if (CreateDirectory(folder, NULL))
	{
		//std::cout << "Directory Created!" << endl;
		return true;
	}
	else if (ERROR_ALREADY_EXISTS == GetLastError())
	{
		//std::cout << "Directory Already Exists!" << endl;
		return true;
	}
	//std::cout << "Error at Creating Directory" << endl;
	return false;

	delete[] folder;
}

#pragma endregion

//void ToolBox::GetDataEvent(double*** seismicConditioning, Coord seismicDataEvent, double*** seismicMatrix, ConfigFile* configFile)
//{
//	int matrixWidth = configFile->GetTemplateWidth();
//	int matrixHeight = configFile->GetTemplateHeight();
//	int matrixDepth = configFile->GetTemplateDepth();
//
//	int x = seismicDataEvent.x;
//	int y = seismicDataEvent.y;
//	int z = seismicDataEvent.z;
//
//	for (int d = 0; d < matrixDepth; d++)
//		for (int h = 0; h < matrixHeight; h++)
//			for (int w = 0; w < matrixWidth; w++)
//				seismicMatrix[w][h][d] = seismicConditioning[x + w][y + h][z + d];
//
//}


//template <class T>
void ToolBox::PrintMatrix(short** matrix, int width, int height)
{
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl << endl;
}

void ToolBox::PrintMatrix(double** matrix, int width, int height)
{
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl << endl;
}






template <class T>
void ToolBox::MEBC_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, T** answer)// , ConfigFile* configFile)
{
	// STEP 1. Define an error surface, having a rectangular size TEMPLATESIZE by OVERLAP.
	double** surfaceError = InitializeMatrix_double(sizeX_template, sizeY_overlap);		// Create the surface error
	InitializeSurfaceError_2D(AMatrix, BMatrix, sizeX_template, sizeY_overlap, surfaceError);

	// STEP 2. Compute the cumulative minimum error along the cutting direction
	double** cumulMinError = InitializeMatrix_double(sizeX_template, sizeY_overlap);// Create the Cumulative Minimum Error table.
	InitializeCumulativeMinimumErrorSurface_2D(surfaceError, sizeX_template, sizeY_overlap, cumulMinError);
	DeleteMatrix(surfaceError, sizeX_template, sizeY_overlap);

	// STEP 3. and STEP 4.
	int * boundary = new int[sizeX_template];// defining the *boundary*
	FindBoundaryHairline(cumulMinError, sizeX_template, sizeY_overlap, boundary);
	DeleteMatrix(cumulMinError, sizeX_template, sizeY_overlap);

	// STEP 5. Merge the two Surfaces (Matrices) using the calculated Boundary.
	MergeMatricesUsingBoundaryCut_2D(AMatrix, BMatrix, sizeX_template, sizeY_overlap, boundary, answer);
	delete[] boundary;

}


template <class T>
void ToolBox::MEBC_3D(T*** AMatrix, T*** BMatrix, int sizeX_template, int sizeY_template, int sizeZ_overlap, T*** answer)
{
	// STEP 1. Define an error surface, having size TEMPLATESIZE x TEMPLATESIZE x OVERLAP.
	double*** surfaceError = InitializeMatrix_double(sizeX_template, sizeY_template, sizeZ_overlap);		// Create the surface error
	InitializeSurfaceError_3D(AMatrix, BMatrix, sizeX_template, sizeY_template, sizeZ_overlap, surfaceError);


	// STEP 2. Compute the cumulative minimum error along the cutting direction:
	double*** cumulMinError = InitializeMatrix_double(sizeX_template, sizeY_template, sizeZ_overlap);	// Create the Cumulative Minimum Error table.
	InitializeCumulativeMinimumErrorSurface_3D(surfaceError, sizeX_template, sizeY_template, sizeZ_overlap, cumulMinError);
	DeleteMatrix(surfaceError, sizeX_template, sizeY_template, sizeZ_overlap);


	// STEP 3. and STEP 4.
	int** boundary = InitializeMatrix_int(sizeX_template, sizeY_template);
	FindBoundarySurface(cumulMinError, sizeX_template, sizeY_template, sizeZ_overlap, boundary);
	DeleteMatrix(cumulMinError, sizeX_template, sizeY_template, sizeZ_overlap);


	// STEP 5. Merge the two Surfaces (Matrices) using the calculated Boundary.
	MergeMatricesUsingBoundaryCut_3D(AMatrix, BMatrix, sizeX_template, sizeY_template, sizeZ_overlap, boundary, answer);
	DeleteMatrix(boundary, sizeX_template, sizeY_template);
}








int ToolBox::IndexOfSmallestValue(double* vector, int vectorSize)
{
	int smallestIndex = 0;						// Here, the index of the smallest value will be stored.

	for (int j = 1; j < vectorSize; j++)		// we iterate over the values...
		if (vector[j] < vector[smallestIndex])	// and if there's a smaller value than the smallest known...
			smallestIndex = j;					// Improvise. Adapt. Overcome.  Or simply, memorize the new smallest.

	return smallestIndex;						// return the index of that tiny dude over there.
}


template<class T>
void ToolBox::InitializeSurfaceError_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, double** surfaceError)
{
	double ans = 0;									// The one that will be multiplied by itself.

	for (int i = 0; i < sizeX_template; i++)		// for each row...
	{
		for (int j = 0; j < sizeY_overlap; j++)		// and each column...
		{
			ans = AMatrix[i][j] - BMatrix[i][j];		// get the value here

														// facieA = AMatrix[i][j];					// leaving this here, for when the WeightMatrix	implementation comes...
														// facieB = BMatrix[i][j];
														// ans = WeightMatrix[facieA][facieB];

			surfaceError[i][j] = ans * ans;			// exquisite, excellent initialization!
		}
	}
}

template<class T>
void ToolBox::InitializeSurfaceError_3D(T*** AMatrix, T*** BMatrix, int sizeX_template, int sizeY_template, int sizeZ_overlap, double*** surfaceError)
{
	double ans = 0;									// The one that will be multiplied by itself.

	for (int i = 0; i < sizeX_template; i++)		// for each layer...
	{
		for (int j = 0; j < sizeY_template; j++)	// and each row...
		{
			for (int k = 0; k < sizeZ_overlap; k++)	// and each column...
			{
				ans = AMatrix[i][j][k] - BMatrix[i][j][k];	// get the value here

															// facieA = AMatrix[i][j][k];				// leaving this here, for when the WeightMatrix	implementation comes...
															// facieB = BMatrix[i][j][k];
															// ans = WeightMatrix[facieA][facieB];

				surfaceError[i][j][k] = ans * ans;	// excellent initialization.
			}
		}
	}
}



void ToolBox::InitializeCumulativeMinimumErrorSurface_2D(double** surfaceError, int sizeX_template, int sizeY_overlap, double** cumulMinError)
{
	double leMinimum = 0;

	for (int i = 0; i < sizeX_template; i++)
	{
		for (int j = 0; j < sizeY_overlap; j++)
		{
			if (i > 0)
			{
				// | B |
				// | A |
				// | C |
				/* A */				leMinimum = cumulMinError[i - 1][j];							// initialize as the value in the previous row, same column
				if (j > 0)
					/* B */				leMinimum = min(leMinimum, cumulMinError[i - 1][j - 1]);	// if there's a anterior column, consider it.
				if (j < sizeY_overlap - 1)
					/* C */				leMinimum = min(leMinimum, cumulMinError[i - 1][j + 1]);	// if there's a posterior column, consider it.

				cumulMinError[i][j] = surfaceError[i][j] + leMinimum;			// aaaaand finally, the Cumulative Minimum Error is the current value + the minimum.
			}
			else
				cumulMinError[i][j] = surfaceError[i][j];
		}
	}
}

void ToolBox::InitializeCumulativeMinimumErrorSurface_3D(double*** surfaceError, int sizeX_template, int sizeY_template, int sizeZ_overlap, double*** cumulMinError)
{
	double leMinimum = 0;

	for (int i = 0; i < sizeX_template; i++)				// for each layer 'i' in 'e' (with 'i' = 2...p), do:
	{
		for (int j = 0; j < sizeY_template; j++)				// for each voxel [j,k] in layer 'i', do:
		{
			for (int k = 0; k < sizeZ_overlap; k++)
			{
				// Calculate the cumulative minimum error E using the 9 closest voxels on the previous layer k-1:
				if (i > 0)
				{
					// From the previous 'layer' (the i-1 one), we find the minimum value.
					//  | E | B | H |
					//	| D | A | G |
					//	| F | C | I |

					/* A */						leMinimum = cumulMinError[i - 1][j][k];						// initialize as the value in the previous layer, same row and column
					if (k > 0)
						/* B */						leMinimum = min(cumulMinError[i - 1][j][k - 1], leMinimum);	// if there's a anterior row, consider it.
					if (k < sizeZ_overlap - 1)
						/* C */						leMinimum = min(cumulMinError[i - 1][j][k + 1], leMinimum);	// if there's a posterior row, consider it.

					if (j > 0)
					{
						/* D */						leMinimum = min(cumulMinError[i - 1][j - 1][k], leMinimum);	// if there's a anterior column, consider it.
						if (k > 0)
							/* E */						leMinimum = min(cumulMinError[i - 1][j - 1][k - 1], leMinimum);	// if there's a anterior row, consider it.
						if (k < sizeZ_overlap - 1)
							/* F */						leMinimum = min(cumulMinError[i - 1][j - 1][k + 1], leMinimum);	// if there's a posterior row, consider it.
					}
					if (j < sizeY_template - 1)
					{
						/* G */						leMinimum = min(cumulMinError[i - 1][j + 1][k], leMinimum);	// if there's a posterior column, consider it.
						if (k > 0)
							/* H */						leMinimum = min(cumulMinError[i - 1][j + 1][k - 1], leMinimum);	// if there's a anterior row, consider it.
						if (k < sizeZ_overlap - 1)
							/* I */						leMinimum = min(cumulMinError[i - 1][j + 1][k + 1], leMinimum);	// if there's a posterior row, consider it.

					}
					cumulMinError[i][j][k] = surfaceError[i][j][k] + leMinimum;			// aaaaand finally, the Cumulative Minimum Error is the current value + the minimum.
				}
				else
					cumulMinError[i][j][k] = surfaceError[i][j][k];
			}
		}
	}
}





void ToolBox::FindBoundaryHairline(double** cumulMinError, int sizeX_template, int sizeY_overlap, int* boundary)
{
	// STEP 3. Identify the coordinate J corresponding to the entry with smallest value on the last row of E. This location	corresponds to the arrival point of a path of minimum cost through the error surface.
	int currentI = sizeX_template - 1;										// the analysis starts in the last row of the Cumulative Minimum Error table.
	int smallestJIndex = IndexOfSmallestValue(cumulMinError[currentI], sizeY_overlap);	// and this guy will hold the column index of the smallest value in the current row.

	boundary[currentI] = smallestJIndex;									// and don't forget to put the smallest guy in the last position of the boundary!!!


																			// STEP 4. Trace back the minimum values for each row i, going backward (with i=p-1...1) and each time identify the cutting path.
	int previousJIndex;														// now, this 'temp' guy will facilitate the readability of the code.

	for (currentI--; currentI >= 0; currentI--)								// iterating over the rest of the rows, in reverse order...
	{
		// From the previous 'layer' (the i-1 one), we find the K-index that provides minimum value.
		// | B |
		// | A |
		// | C |

		previousJIndex = boundary[currentI + 1];							// get the column index of the   previous row smallest value
		/* A */	smallestJIndex = boundary[currentI + 1];							// and consider it as a possible current  row smallest value

																					// if the smallest of the previous row had a column to the 'left' of it, and if the value in this 'left' column is smaller than the smallest known...
		if (previousJIndex > 0 && cumulMinError[currentI][previousJIndex - 1] < cumulMinError[currentI][smallestJIndex])
			/* B */		smallestJIndex = previousJIndex - 1;

		// if the smallest of the previous row had a column to the 'right' of it, and if the value in this 'right' column is smaller than the smallest known...
		if (previousJIndex < sizeY_overlap - 1 && cumulMinError[currentI][previousJIndex + 1] < cumulMinError[currentI][smallestJIndex])
			/* C */		smallestJIndex = previousJIndex + 1;

		boundary[currentI] = smallestJIndex;								// Include the column index of the smallest value of this row in the boundary.
	}
}

void ToolBox::FindBoundarySurface(double*** cumulMinError, int sizeX_template, int sizeY_template, int sizeZ_overlap, int** boundary)
{
	// STEP 3. Consider only the last layer of E. Each E value on this last layer represents the cost of the shortest 1-D
	// path throughout the block. We isolate this last layer and perform a 2-D cut (section 2.2) to find the minimum
	// error cut in E_(p-1),i,j. This 2-D cut encompasses the maximum number of arrival points of 1-D shortest
	// paths throughout the block.
	int currentI = sizeX_template - 1;										// the analysis starts in the last row of the Cumulative Minimum Error table.
	FindBoundaryHairline(cumulMinError[currentI], sizeY_template, sizeZ_overlap, boundary[currentI]);


	int previousK_index;
	int smallestK_index;


	// STEP 4. Propagate this 2-D cut throughout the block by going backward (with k=p-1...1)
	for (currentI--; currentI >= 0; currentI--)								// iterating over the rest of the rows, in reverse order...
	{
		for (int currentJ = 0; currentJ < sizeY_template; currentJ++)				// for each voxel [j,k] in layer 'i', do:
		{
			// From the previous 'layer' (the i-1 one), we find the K-index that provides minimum value.
			// | B |
			// | A |
			// | C |

			previousK_index = boundary[currentI + 1][currentJ];							// get the column index of the   previous row smallest value
			/* A */		smallestK_index = boundary[currentI + 1][currentJ];							// and consider it as a possible current  row smallest value

																									// if the smallest of the previous row had a column to the 'left' of it, and if the value in this 'left' column is smaller than the smallest known...
			if (previousK_index > 0 && cumulMinError[currentI][currentJ][previousK_index - 1] < cumulMinError[currentI][currentJ][smallestK_index])
				/* B */			smallestK_index = previousK_index - 1;

			// if the smallest of the previous row had a column to the 'right' of it, and if the value in this 'right' column is smaller than the smallest known...
			if (previousK_index < sizeZ_overlap - 1 && cumulMinError[currentI][currentJ][previousK_index + 1] < cumulMinError[currentI][currentJ][smallestK_index])
				/* C */			smallestK_index = previousK_index + 1;

			boundary[currentI][currentJ] = smallestK_index;								// Include the column index of the smallest value of this row in the boundary.
		}
	}
}





template<class T>
void ToolBox::MergeMatricesUsingBoundaryCut_2D(T** AMatrix, T** BMatrix, int sizeX_template, int sizeY_overlap, int* boundary, T** answer)
{
	for (int i = 0; i < sizeX_template; i++)
	{
		for (int j = 0; j < boundary[i]; j++)
			answer[i][j] = AMatrix[i][j];
		for (int j = boundary[i]; j < sizeY_overlap; j++)
			answer[i][j] = BMatrix[i][j];
	}
}

template<class T>
void ToolBox::MergeMatricesUsingBoundaryCut_3D(T*** AMatrix, T*** BMatrix, int sizeX_template, int sizeY_template, int sizeZ_overlap, int** boundary, T*** answer)
{
	for (int i = 0; i < sizeX_template; i++)
	{
		for (int j = 0; j < sizeY_template; j++)
		{
			for (int k = 0; k < boundary[i][j]; k++)
				answer[i][j][k] = AMatrix[i][j][k];
			for (int k = boundary[i][j]; k < sizeZ_overlap; k++)
				answer[i][j][k] = BMatrix[i][j][k];
		}
	}
}