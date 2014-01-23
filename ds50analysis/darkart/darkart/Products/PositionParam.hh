/*
 Definiton of the PositionParam class

 Holds information about the 3D position reconstruction of
 an event. 
*/

#ifndef darkart_Products_PositionParam_hh
#define darkart_Products_PositionParam_hh

#include <vector>

namespace darkart {
  struct PositionParam;
  typedef std::vector<PositionParam> PositionParamVec;
}


struct darkart::PositionParam
{
  double bary_x;
  double bary_y;

  PositionParam():
    bary_x(0.),
    bary_y(0.)
  { }

};

inline bool operator == (darkart::PositionParam const& a, darkart::PositionParam const& b)
{
  return (a.bary_x == b.bary_x &&
          a.bary_y == b.bary_y );
}




#endif
