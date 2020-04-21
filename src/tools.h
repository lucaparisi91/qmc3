#ifndef TOOLS_H
#define TOOLS_H

#include "traits.h"
#include <fstream>

std::string ansiColor(const std::string & color_name);

bool is_empty(std::ifstream& pFile);


inline int getN(const state_t & state){  return state.rows();};

inline size_t size(const distance_t & dis) {return dis.rows();}

inline constexpr int getDimensions() { return DIMENSIONS; }





//inline constexpr int getDimensions(const state_t & state) {return DIMENSIONS;}
#endif
