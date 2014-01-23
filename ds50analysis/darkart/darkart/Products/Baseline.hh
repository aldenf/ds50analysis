/*

  Definition of the Baseline class

  AFan 2013-07-03

 */


#ifndef darkart_Products_Baseline_hh
#define darkart_Products_Baseline_hh

#include <vector>
#include <map>
#include <utility>
#include <iostream>

namespace darkart
{
  // Forward declarations
  struct Baseline;
  typedef std::vector<Baseline> BaselineVec;
}

  
// Holds useful data about the baseline for a single
// channel as found by the baseline finder
struct darkart::Baseline
{
  // Default constructor
  Baseline():
    found_baseline(false),
    mean(-1.0),
    variance(-1.0),
    saturated(false),
    laserskip(false)
  { }

  // Useful info about the baseline
  bool    found_baseline; //Was the baseline finder successful?
  double            mean; //the average baseline found
  double        variance; //variance of the samples in the baseline range

  bool         saturated;
  int             length;
  int search_start_index;

  std::vector< std::pair< double, double > > interpolations; //start and end time of each baseline interpolation
  
  bool         laserskip;

  bool isValid() const {
    return (found_baseline && mean>-1 && variance>-1 &&
            length>-1 && search_start_index>-1);
  }

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "Baseline Info:"                                              << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << (found_baseline  ? "Baseline found." : "Baseline not found.") << std::endl;
    std::cout << "Mean:                       " << mean                        << std::endl;
    std::cout << "Variance:                   " << variance                    << std::endl;

    if(verbosity > 0){
      std::cout << "Saturated:                  " << saturated                   << std::endl;
      std::cout << "Length:                     " << length                      << std::endl;
      std::cout << "Search Start Index:         " << search_start_index          << std::endl;
      std::cout << "Laserskip:                  " << laserskip                   << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }

};

inline
bool operator==(darkart::Baseline const& a, darkart::Baseline const& b)
{
  return
    a.found_baseline == b.found_baseline &&
    a.mean == b.mean &&
    a.variance == b.variance &&
    a.saturated == b.saturated;
}

#endif
