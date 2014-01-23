/*

  This file is necessary to create the dictionary used to build
  the ROOT tree. This is necessary for every class that is an
  instantiation of a class template (e.g., std::vector<darkart::BaselineData>
  is an instantiation of the class template std::vector). Note
  this is necessary for the art::Wrapper instantiation for each
  top-level data product.

  For most data products, you may mimic the pattern of the
  BaselineData object below:
  1. wrapper for the object
  2. vector of objects
  3. wrapper for the vector of objects
  4. (as necessary) wrapper for the Assns product

  For many cases, the wrapper for the bare object (1) is not
  necessary; the vector of objects will be used. But add an entry
  for (1) anyway to avoid possible nasty debugging situations.

  Make sure to #include the corresponding header file.

 */


#include "darkart/Products/EventData.hh"
#include "darkart/Products/EventInfo.hh"
#include "darkart/Products/Waveform.hh"
#include "darkart/Products/WaveformInfo.hh"
#include "darkart/Products/Channel.hh"
#include "darkart/Products/Baseline.hh"
#include "darkart/Products/Pmt.hh"
#include "darkart/Products/ProductAssns.hh"
#include "darkart/Products/VetoTruth.hh"
#include "darkart/Products/LaserSpectrum.hh"
#include "darkart/Products/Roi.hh"
#include "darkart/Products/Pulse.hh"
#include "darkart/Products/PulseParam.hh"
#include "darkart/Products/PositionParam.hh"
#include "darkart/Products/Spe.hh"
#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"
#include "art/Persistency/Common/Ptr.h"
#include <vector>
#include <utility>

template class std::vector< std::map< double, double> >;

template class art::Wrapper<darkart::EventInfo>;

template class art::Wrapper<darkart::Channel>;
template class std::vector<darkart::Channel>;
template class art::Wrapper<std::vector<darkart::Channel> >;

template class art::Wrapper<darkart::Channel::ChannelID>;

template class art::Wrapper<darkart::Waveform>;
template class std::vector<darkart::Waveform>;
template class art::Wrapper<std::vector<darkart::Waveform> >;
template class art::Wrapper<darkart::WaveformAssns>;

template class art::Wrapper<darkart::WaveformInfo>;
template class std::vector<darkart::WaveformInfo>;
template class art::Wrapper<std::vector<darkart::WaveformInfo> >;
template class art::Wrapper<darkart::WaveformInfoAssns>;

template class art::Wrapper<darkart::Pmt>;
template class std::vector<darkart::Pmt>;
template class art::Wrapper<std::vector<darkart::Pmt> >;

template class art::Wrapper<darkart::Baseline>;
template class std::vector<darkart::Baseline>;
template class art::Wrapper<std::vector<darkart::Baseline> >;
template class art::Wrapper<darkart::BaselineAssns>;

template class art::Wrapper<darkart::Roi>;
template class std::vector<darkart::Roi>;
template class art::Wrapper<std::vector<darkart::Roi> >;
template class art::Wrapper<darkart::RoiAssns>;

template class art::Wrapper<darkart::Pulse>;
template class art::Wrapper<std::vector<darkart::Pulse> >;
template class art::Wrapper<darkart::PulseAssns>;
template class std::vector<darkart::Pulse>;

template class art::Wrapper<darkart::Pulse::PulseID>;

template class art::Wrapper<darkart::PulseParam>;
template class art::Wrapper<std::vector<darkart::PulseParam> >;
template class art::Wrapper<darkart::PulseParamAssns>;
template class std::vector<darkart::PulseParam>;

template class art::Wrapper<darkart::PositionParam>;
template class art::Wrapper<std::vector<darkart::PositionParam> >;
template class art::Wrapper<darkart::PositionParamAssns>;
template class std::vector<darkart::PositionParam>;

template class art::Wrapper<darkart::Spe>;
template class art::Wrapper<std::vector<darkart::Spe> >;
template class art::Wrapper<darkart::SpeAssns>;
template class std::vector<darkart::Spe>;

template class std::vector<darkart::VetoTDCHit>;
template class art::Wrapper<darkart::VetoTruth>;

template class art::Wrapper<darkart::LaserSpectrum>;
template class art::Wrapper<TFitResultPtr>;
template class art::Wrapper<darkart::LaserResults>;
template class cet::map_vector<darkart::LaserSpectrum*>;
template class art::Wrapper<cet::map_vector<darkart::LaserSpectrum*> >;
template class std::vector<std::pair<cet::map_vector_key,darkart::LaserSpectrum*> >;
template class art::Wrapper<std::vector<std::pair<cet::map_vector_key,darkart::LaserSpectrum*> > >;
template class std::pair<cet::map_vector_key,darkart::LaserSpectrum*>;
template class art::Wrapper<std::pair<cet::map_vector_key,darkart::LaserSpectrum*> >;



template class std::vector<darkart::PulseData>;
template class std::vector<darkart::ChannelData>;
template class art::Wrapper<darkart::EventData>;
