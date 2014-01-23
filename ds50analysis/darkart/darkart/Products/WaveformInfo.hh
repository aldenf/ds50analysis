/*

  Definition of the WaveformInfo class

  AFan 2013-07-29

 */


#ifndef darkart_Products_WaveformInfo_hh
#define darkart_Products_WaveformInfo_hh

#include <vector>
#include <iostream>

// Forward declarations
namespace darkart
{
  struct WaveformInfo;
  typedef std::vector<WaveformInfo> WaveformInfoVec;
}


struct darkart::WaveformInfo
{
  WaveformInfo():
    minimum(0.),
    maximum(0.),
    min_index(-1),
    max_index(-1),
    min_time(-1.),
    max_time(-1.),
    nsamps(-1)
  { }

  double       minimum; //minimum value of a waveform
  double       maximum; //maximum value of a waveform
  int        min_index; //index at which minimum value occurs in waveform
  int        max_index; //index at which maximum value occurs in waveform
  double      min_time; //time at which minimum value occurs in waveform
  double      max_time; //time at which maximum value occurs in waveform
  int           nsamps; //number of samples in the waveform

  bool isValid() const {
    return (min_index>-1 && max_index>-1);
  }

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "Waveform Info (details):"                                    << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << "Minimum:                    " << minimum                     << std::endl;
    std::cout << "Maximum:                    " << maximum                     << std::endl;
    std::cout << "Number of Samples:          " << nsamps                      << std::endl;

    if(verbosity > 0){
      std::cout << "Index of Minimum:           " << min_index                   << std::endl;
      std::cout << "Time of Minimum:            " << min_time                    << std::endl;
      std::cout << "Index of Maximum:           " << max_index                   << std::endl;
      std::cout << "Time of Maximum:            " << max_time                    << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }
  
};

inline
bool operator==(darkart::WaveformInfo const& a, darkart::WaveformInfo const& b)
{
  return
    a.minimum == b.minimum &&
    a.maximum == b.maximum &&
    a.min_index == b.min_index &&
    a.max_index == b.max_index &&
    a.min_time == b.min_time &&
    a.max_time == b.max_time &&
    a.nsamps == b.nsamps;
}



#endif
