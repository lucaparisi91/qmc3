#include "tools.h"
#include "qmcExceptions.h"

std::string ansiColor(const std::string & color_name)
{
	if (color_name == "red") return  "\033[0;31m"; 
	if (color_name == "default") return  "\033[0m";
	if (color_name == "cyan") return  "\033[0;36m"; 
	if (color_name == "green") return  "\033[0;32m"; 

	throw invalidInput(color_name);

};

bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}
