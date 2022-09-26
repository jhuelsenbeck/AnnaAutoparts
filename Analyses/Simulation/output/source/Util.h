#ifndef Util_H
#define Util_H

#include <iostream>
#include <string>
#include <fstream>



std::string getLineFromFile(std::string fileName, int lineNum) {

	/* open the file */
	std::ifstream fileStream(fileName.c_str());
	if (!fileStream) 
		{
		std::cerr << "Cannot open file \"" + fileName + "\"" << std::endl;
		exit(1);
		}

	std::string linestring = "";
	int line = 0;
	while( getline(fileStream, linestring).good() )
		{
		line++;
		if (line == lineNum)
			break;
		}

	/* close the file */
	fileStream.close();
	
	if (line != lineNum)
		{
		std::cerr << "The file \"" + fileName + "\" has " << line << " lines. Could not find line " << lineNum << std::endl;
		exit(1);
		}
	
	return linestring;
}

int flip(int x) {
	
	if (x == 0)
		return 1;
	else
		return 0;
		
}

#endif