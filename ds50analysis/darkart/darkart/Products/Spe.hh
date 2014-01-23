/* 
Definition of the Spe class
Masayuki Wada 2013-11-08
*/

#ifndef darkart_Products_Spe_hh
#define darkart_Products_Spe_hh

#include <vector>
#include <iostream>

namespace darkart
{
  struct Spe;
  typedef std::vector<Spe> SpeVec;
}

struct darkart::Spe
{
  Spe():
    integral(0.),
    start_time(0.),
    amplitude(0.),
    peak_time(0.),
    local_baseline(0.),
    length(0.),
    start_clean(true),
    end_clean(true),
    dt(0.)
  {}
  double integral; ///< Integral in counts*samples of the found charge
  double start_time; ///< time of the start of the pulse
  double amplitude;  ///< Maximum height obtained relative to local baseline
  double peak_time; ///< time at which the pulse reached max amplitude
  double local_baseline; ///< value of the local baseline before the pulse
  double length; ///< Time in us over which integral was evaluated
  bool start_clean;    ///< start of spe does not overlap with previous
  bool end_clean;      ///< end of spe does not overlap with next
  double dt;     ///< Time in us since the previous found spe

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "Single Photoelectron (SPE) Info:"                            << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << "Integral:                   " << integral                    << std::endl;
    std::cout << "Amplitude:                  " << amplitude                   << std::endl;
    std::cout << "Time of Peak:               " << peak_time                   << std::endl;

    if(verbosity > 0){
      std::cout << "Pulse Start Time:         " << start_time                  << std::endl;
      std::cout << "Duration of Pulse:        " << length                      << std::endl;
      std::cout << "Local Baseline:           " << local_baseline              << std::endl;
      std::cout << "Clean Start:              " << start_clean                 << std::endl;
      std::cout << "Clean End:                " << end_clean                   << std::endl;
      std::cout << "Time Since Previous:      " << dt                          << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }
    
};

inline bool operator== (darkart::Spe const& a, darkart::Spe const& b)
{
  return (a.integral == b.integral
      && a.start_time == b.start_time
      && a.amplitude == b.amplitude && a.peak_time == b.peak_time
      && a.local_baseline == b.local_baseline && a.length == b.length
      && a.start_clean == b.start_clean && a.end_clean == b.end_clean
      && a.dt == b.dt);
}

#endif
