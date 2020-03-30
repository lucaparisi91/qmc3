#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include "traits.h"

struct parameter
{
  parameter( std::string name_="", size_t  offset_=0,size_t length_=0) : name(name_),offset(offset_),length(length_) {};

  template<class T>
  auto begin(T & storage){return storage.begin() + offset;};

  template<class T>
  auto cbegin(T & storage) const {return storage.cbegin() + offset;};

  template<class T>
  auto end(T & storage){return storage.begin() + offset + length;}

  template<class T>
  auto end(T & storage) const {return storage.cbegin() + offset + length;}

  auto size() const {return length;}

  size_t offset;
  size_t length;
  std::string name;
};

#endif