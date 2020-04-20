#include "abstractFactory.h"

std::string createId(const json_t & j)
  {
    std::string id="";
    if ( j.find("kind") != j.end() )
      id +=j["kind"];
    
    for (auto  & element : j)
      {
	if ( element.is_structured() )
	  {
	    id+="/"+createId(element);
	  }

      }
    return id;
  }

