#include "FileReadingToolBox.h"



string FileReadingToolBox::GetNextValidLine(ifstream& fs)
{
	for (string line; getline(fs, line); )
	{
		// TODO: TRIM THE LINE'S BEGINNING AND ENDING!!!!!!!!!!!!!!!!!!!!!!!


		if (line[0] == '#')
			continue;

		if (line.empty())
			continue;

		return line;
	}
	return "";
}


vector<string> FileReadingToolBox::SplitString(const string& parameterString, const string& delimiter)
{
	string myString = parameterString;
	vector<string> tokens;

	size_t pos = 0;
	string token;
	while ((pos = myString.find(delimiter)) != string::npos) {
		token = myString.substr(0, pos);
		tokens.push_back(token);
		myString.erase(0, pos + delimiter.length());
	}
	tokens.push_back(myString);
	return tokens;
}
