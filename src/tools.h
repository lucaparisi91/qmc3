#ifndef TOOLS_H
#define TOOLS_H

#include "traits.h"
#include <fstream>
#include <vector>


json_t toJson(const states_t & states);
std::vector<states_t> readStates(json_t & jI);

std::string ansiColor(const std::string & color_name);

bool is_empty(std::ifstream& pFile);


inline int getN(const state_t & state){  return state.rows();};

inline size_t size(const distance_t & dis) {return dis.rows();}

inline constexpr int getDimensions() { return DIMENSIONS; }

std::vector<states_t> readStatesFromDirectory(std::string filename);



int wrapIndex(int i , int size);

//inline constexpr int getDimensions(const state_t & state) {return DIMENSIONS;}
#endif
