#ifndef darkart_Products_VetoTruth_hh
#define darkart_Products_VetoTruth_hh

#include <vector>

namespace darkart
{
	// An instance of VetoTDCHit represents the MC truth information for
	// a single hit.
  struct VetoTDCHit
  {
    int   pmt_index;
    int   hit_index;
    float leading_ns;
    float trailing_ns;
    float pmt_theta;
    float pmt_phi;
  };

  typedef std::vector<VetoTDCHit> VetoTDCHits;

  struct VetoTruth
  {
    int         sim_event;
		VetoTDCHits hits;
  };

}

#endif

