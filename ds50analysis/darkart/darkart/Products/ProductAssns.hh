#ifndef darkart_Products_ProductHierarchy_h
#define darkart_Products_ProductHierarchy_h

/*
  Relatioships between data products. Add an entry for
  each channel-level data product.
 */

#include "art/Persistency/Common/Assns.h"

#include "darkart/Products/Channel.hh"
#include "darkart/Products/Baseline.hh"
#include "darkart/Products/Waveform.hh"
#include "darkart/Products/WaveformInfo.hh"
#include "darkart/Products/Roi.hh"
#include "darkart/Products/Pulse.hh"
#include "darkart/Products/PulseParam.hh"
#include "darkart/Products/PositionParam.hh"
#include "darkart/Products/Spe.hh"

namespace darkart {
  typedef art::Assns<Channel, Baseline> BaselineAssns;
  typedef art::Assns<Channel, Waveform> WaveformAssns;
  typedef art::Assns<Channel, WaveformInfo> WaveformInfoAssns;
  typedef art::Assns<Channel, Roi> RoiAssns;
  typedef art::Assns<Channel, Pulse> PulseAssns;
  typedef art::Assns<Pulse, PulseParam> PulseParamAssns;
  typedef art::Assns<Pulse, PositionParam> PositionParamAssns;
  typedef art::Assns<Channel, Spe> SpeAssns;
}


#endif
