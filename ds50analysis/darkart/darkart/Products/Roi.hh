/*
Definition of the RoiData class
AKurlej 2013-08-01
 */

#ifndef darkart_Products_RoiData_hh
#define darkart_Products_RoiData_hh

#include <vector>
#include <iostream>

namespace darkart
{
  struct Roi;
  typedef std::vector<Roi> RoiVec;

}


struct darkart::Roi
{
  //constructors n' functions, yo
  Roi(): 
    region_number(-1), start_time(0), end_time(0),
    start_index(-1), end_index(-1),
    max(0),max_time(0), min(0), min_time(0),
    integral(0)
  { }

  //info
  int region_number;
  double start_time;//input parameter
  double end_time;//input parameter
  int start_index;
  int end_index;
  double max;
  double max_time;
  double min;
  double min_time;  
  double integral;
  //double npe;
  //int min_index;

  bool isValid() const {
    return (region_number>-1 && start_index>-1 && end_index>-1);
  }

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "Region of Interest (ROI) Info:"                              << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << "Region Number:              " << region_number               << std::endl;
    std::cout << "Start Time:                 " << start_time                  << std::endl;
    std::cout << "End Time:                   " << end_time                    << std::endl;

    if(verbosity > 0){
      std::cout << "Max over ROI:               " << max                         << std::endl;
      std::cout << "Time of Max:                " << max_time                    << std::endl;
      std::cout << "Min over ROI:               " << min                         << std::endl;
      std::cout << "Time of Min:                " << min_time                    << std::endl;
      std::cout << "Integral:                   " << integral                    << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }
  
};

inline
bool operator==(darkart::Roi const& a, darkart::Roi const& b)
{
return
  a.start_time == b.start_time &&
  a.end_time == b.end_time &&
  a.max == b.max &&
  a.max_time == b.max_time &&
  a.min == b.min &&
  a.min_time == b.min_time &&
  a.integral == b.integral;
}

#endif
