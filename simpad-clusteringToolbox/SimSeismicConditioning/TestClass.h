#ifndef TEST_CLASS
#define TEST_CLASS

#include <iostream>
#include "ToolBox.h"

using namespace std;

class TestClass
{
public:
	static void TestAll();

	static void TestStringSplit();


};


void TestClass::TestAll()
{
	cout << "TestClass::TestAll()" << endl;
	TestStringSplit();

}

void TestClass::TestStringSplit()
{
	cout << "Beginning TestStringSplit:" << endl;
	bool testFailed = false;

	string myString = "uno dos tres cuatro";
	vector<string> originalTokens = { "uno", "dos", "tres", "cuatro" };
	string delimiter = " ";

	vector<string> tokens = FileReadingToolBox::SplitString(myString, delimiter);


	if (originalTokens.size() == tokens.size())
	{
		for (int i = 0; i < originalTokens.size(); i++)
		{
			if (originalTokens[i] != tokens[i])
			{
				cout << "Test Failed: Different tokens:" << originalTokens[i] << "," << tokens[i] << endl;
				testFailed = true;
			}
		}
	}
	else
	{
		cout << "Test failed: Incorrect amount of tokens" << endl;
		testFailed = true;
	}

	if (!testFailed)
		cout << "Test TestStringSplit SUCCEEDED" << endl;
}




#endif // TEST_CLASS
