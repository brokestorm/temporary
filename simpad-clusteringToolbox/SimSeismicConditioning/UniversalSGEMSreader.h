#ifndef UNIVERSAL_SGEMS_READER
#define UNIVERSAL_SGEMS_READER

#include <string>
//#include <iostream>
#include <fstream>
#include <vector>


#include "FileReadingToolBox.h"

using namespace std;

class UniversalSGEMSreader
{
public:
	static short*** ReadSGEMS_3D_short1(string path, int& auxW, int& auxH, int& auxD);
	static short*** ReadSGEMS_3D_short_FAST(string path, int& auxW, int& auxH, int& auxD);
	static double*** ReadSGEMS_3D_double_FAST(string path, int& auxW, int& auxH, int& auxD);

private:
	static short*** InitializeMatrix_short(int width, int height, int depth);
	static double*** InitializeMatrix_double(int width, int height, int depth);
};

#endif // !UNIVERSAL_SGEMS_READER
