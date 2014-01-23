/* 
Definition of the Pulse class
RSaldanha 2013-09-08
*/

#ifndef darkart_Products_Pulse_hh
#define darkart_Products_Pulse_hh

#include <vector>
#include <iostream>

namespace darkart
{
  struct Pulse;
  typedef std::vector<Pulse> PulseVec;
}

struct darkart::Pulse
{
    struct PulseID 
    {
	int pulse_id;     ///< rank of this pulse within its channel; 1st pulse? 2nd pulse? etc.
	int pulse_index;  ///< position of this pulse within the data product vector
    };
    
    
    PulseID pulseID;     ///< unique identifier for this pulse
    bool start_clean;    ///< start of pulse does not overlap with previous pulse
    bool end_clean;      ///< end of pulse does not overlap with next pulse
    int start_index;     ///< index marking the start of the pulse
    int end_index;       ///< index for the end of the pulse
    double start_time;   ///< time since trigger at start of pulse
    double end_time;     ///< time since trigger at end of pulse
    double dt;           ///< time between the start of this pulse and the previous one (or start of trigger)

    bool is_s1_like; ///< pulse looks like an S1 for purposes of possible start time adjustment
    bool start_shifted;  ///< start of this pulse is shifted or not
    double ratio1;       ///< the fraction of the total pulse integral within +-ratio_samps(20 ns) at peak
    double ratio2;       ///< the fraction of the total pulse integral between the start_index and peak_index-ratio_samps
    double start_shift;  //How much the start was shifted, default = 0 in us

    int pulse_id() const { return pulseID.pulse_id; }
    
    Pulse():
      pulseID(),
      start_clean(false), end_clean(false),
      start_index(-1), end_index (-1),
      start_time(0), end_time(0),
      dt (-1),
      is_s1_like(false),
      start_shifted(false),
      ratio1(0.), ratio2(0.),
      start_shift(0.)
	{}

    // AWW 24.11.13
    inline void print(int verbosity = 0) const {
      std::cout << "Pulse Info:"                                                 << std::endl;
      std::cout << "----------------------------------------"                    << std::endl;
      std::cout << "Pulse ID:           " << pulse_id()                  << std::endl;
      std::cout << "Start Time:         " << start_time                  << std::endl;
      std::cout << "End Time:           " << end_time                    << std::endl;

      if(verbosity > 0){
        std::cout << "Clean Start:        " << start_clean                 << std::endl;
       	std::cout << "Clean End:          " << end_clean                   << std::endl; }
      if(verbosity > 1){
        /* placeholder */ }

      std::cout << std::endl;
    } 

};

inline bool operator== (darkart::Pulse const& a, darkart::Pulse const& b)
{
    return (a.pulseID.pulse_id == b.pulseID.pulse_id &&
	    a.pulseID.pulse_index == b.pulseID.pulse_index &&
	    a.start_clean == b.start_clean &&
	    a.end_clean == b.end_clean &&
	    a.start_index == b.start_index &&
	    a.end_index == b.end_index &&
	    a.start_time == b.start_time &&
	    a.end_time == b.end_time &&
	    a.dt == b.dt
	);
}

#endif
