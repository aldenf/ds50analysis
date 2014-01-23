/*

  Definition of the Waveform class

  AFan 2013-06-30

 */

#ifndef darkart_Products_Waveform_hh
#define darkart_Products_Waveform_hh

#include <vector>
#include <iostream>

namespace darkart
{
  // Forward declarations
  struct Waveform;
  typedef std::vector<Waveform> WaveformVec;
}

  
// A Waveform object holds a full waveform (e.g. a raw waveform or
// baseline-subtracted waveform) for a single channel. 
struct darkart::Waveform
{
  Waveform():
    wave(),
    trigger_index(-1),
    sample_rate(-1.0)
  { }

  std::vector<double> wave;
  int trigger_index;
  double sample_rate; //in MHz

  bool isValid() const;
    
  double SampleToTime(int sample) const;
  int    TimeToSample(double time, bool checkrange=false) const;

  bool operator==(const Waveform & other) const {
    return 
      ( trigger_index == other.trigger_index &&
        sample_rate == other.sample_rate &&
        wave == other.wave );
  }

  bool operator!=(const Waveform & other) const {
    return !(*this == other);
  }

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "Waveform Info:"                                              << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << "Trigger Index:              " << trigger_index               << std::endl;
    std::cout << "Sample Rate (MHz):          " << sample_rate                 << std::endl;

    if(verbosity > 0){ }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }

};


#endif
