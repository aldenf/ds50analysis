/*
Definition of the PulseParam class
*/

#ifndef darkart_Products_PulseParam_hh
#define darkart_Products_PulseParam_hh

#include <map>
#include <utility>
#include <iostream>

namespace darkart{
  struct PulseParam;
  typedef std::vector<PulseParam> PulseParamVec;
}

struct darkart::PulseParam
{
  bool found_peak;       ///< did we find a peak?
  int peak_index;        ///< index marking the peak of the pulse
  double peak_time;      ///< time at which the peak occurs
  double peak_amplitude; ///< amplitude of the peak
  double integral;       ///< integral of the pulse
  bool peak_saturated; ///< was the peak saturated?
  std::vector <double> f_param; ///< f-parameters for different time values
  std::vector <double> f_param_fixed; ///< f-parameters for different time values using fixed_int1
  double f90; ///< f-parameter for 90 ns
  double f90_fixed; ///< f-parameter for 90 ns using fixed_int1
  double t05; ///< time to reach 5% of total integral
  double t10; ///< time to reach 10% of total integral
  double t90; ///< time to reach 90% of total integral
  double t95; ///< time to reach 95% of total integral
  double fixed_int1;  ///< integral evaluated at fixed_time1
  double fixed_int2;  ///< integral evaluated at fixed_time2
  bool fixed_int1_valid; ///< did the event extend at least past fixed_time1?
  bool fixed_int2_valid; ///< did the event extend at least past fixed_time2?
  double npe;            ///< integral scaled for single pe amplitude
  double npe_fixed;      ///< fixed_int1 scaled for single pe amplitude
  double max_ch;         ///< channel ID of the PMT that saw the most light
  std::vector<std::pair<double, double> > baseline_interpolations; ///< start and end times of baseline interpolations in this pulse

  PulseParam():
    found_peak(false), peak_index(-1),
    peak_time(0.), peak_amplitude(0.),
    integral(0.), peak_saturated(false),
    f_param(10), f_param_fixed(10),
    f90(0.), f90_fixed(0.),
    t05(0.), t10(0.),
    t90(0.), t95(0.),
    fixed_int1(0.), fixed_int2(0.),
    fixed_int1_valid(false), fixed_int2_valid(false),
    npe(0.), npe_fixed(0.),
    max_ch(-1),
    baseline_interpolations()
  {}

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "PulseParam Info:"                                            << std::endl;
    std::cout << "Integral:           " << integral                            << std::endl;
    std::cout << "S1 Fixed Integral:  " << fixed_int1                          << std::endl;
    std::cout << "S2 Fixed Integral:  " << fixed_int2                          << std::endl;
    std::cout << "f90:                " << f90                                 << std::endl;

    if(verbosity > 0){
      std::cout << (found_peak  ? "Peak found." : "Peak not found.")             << std::endl;
      std::cout << "Scaled (pe) Integral:       " << npe                         << std::endl;
      std::cout << "Peak Index:                 " << peak_index                  << std::endl;
      std::cout << "Peak Time:                  " << peak_time                   << std::endl;
      std::cout << "Peak Amplitude:             " << peak_amplitude              << std::endl;
      std::cout << (peak_saturated  ? "Peak saturated." : "Peak not saturated.") << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }


};

inline bool operator == (darkart::PulseParam const& a, darkart::PulseParam const& b)
{
    return(a.found_peak == b.found_peak &&
	   a.peak_index == b.peak_index &&
	   a.peak_time == b.peak_time &&
	   a.peak_amplitude == b.peak_amplitude &&
	   a.integral == b.integral &&
	   a.f_param == b.f_param &&
           a.f_param_fixed == b.f_param_fixed &&
	   a.f90 == b.f90 &&
           a.f90_fixed == b.f90_fixed &&
	   a.t05 == b.t05 &&
	   a.t10 == b.t10 &&
	   a.t90 == b.t90 &&
	   a.t95 == b.t95 &&
	   a.fixed_int1 == b.fixed_int1 &&
	   a.fixed_int2 == b.fixed_int2 &&
	   a.fixed_int1_valid == b.fixed_int1_valid &&
	   a.fixed_int2_valid == b.fixed_int2_valid &&
           a.npe == b.npe &&
           a.npe_fixed == b.npe_fixed &&
           a.max_ch == b.max_ch &&
	   a.peak_saturated == b.peak_saturated);
}
#endif
