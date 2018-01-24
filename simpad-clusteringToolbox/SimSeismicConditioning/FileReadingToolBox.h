#ifndef FILE_READING_TOOL_BOX
#define FILE_READING_TOOL_BOX

//#include <cstring>
//#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class FileReadingToolBox
{
public:

	static string FileReadingToolBox::GetNextValidLine(ifstream& fs);

	static vector<string> SplitString(const string& parameterString, const string& delimiter);


private:

};



#endif // ! FILE_READING_TOOL_BOX

