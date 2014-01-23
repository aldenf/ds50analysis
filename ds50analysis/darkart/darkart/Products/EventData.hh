/*

  Defines the EventData data product. Will be a 'super-product'
  which contains all the information from other products. This
  product will be written to a separate output ROOT TTree in a
  separate output file. High-level analysis is done primarily
  on this object.

  AFan 2013-08-23

 */

#ifndef darkart_Products_EventData_hh
#define darkart_Products_EventData_hh

#include "EventInfo.hh"
#include "Channel.hh"
#include "Pmt.hh"
#include "WaveformInfo.hh"
#include "Baseline.hh"
#include "Roi.hh"
#include "Pulse.hh"
#include "PulseParam.hh"
#include "PositionParam.hh"
#include "Spe.hh"

#include <vector>
#include <iostream>
#include <stdexcept>

namespace darkart {
  struct EventData;
  struct ChannelData;
  struct PulseData;
}

struct darkart::PulseData
{
  darkart::Pulse pulse;
  darkart::PulseParam param;
  darkart::PositionParam position;
};

struct darkart::ChannelData
{
  darkart::Channel channel;
  darkart::Pmt pmt;
  darkart::WaveformInfo raw_wfm;
  darkart::WaveformInfo baseline_subtracted_wfm;
  darkart::Baseline baseline;
  std::vector<darkart::Roi> regions;
  darkart::WaveformInfo integral;
  std::vector<darkart::PulseData> pulses;
  std::vector<darkart::Spe> single_pe;
};

struct darkart::EventData
{
  darkart::EventInfo event_info;
  std::vector<darkart::ChannelData> channels;
  darkart::ChannelData sumchannel;

  std::vector<darkart::PulseData> pulses;
  
  void clear();

  darkart::ChannelData const& getChannelByID(int ch_id);
};

inline void darkart::EventData::clear()
{
  event_info = EventInfo();
  channels.clear();
  sumchannel = ChannelData();
  pulses.clear();
}


inline darkart::ChannelData const& darkart::EventData::getChannelByID(int ch_ID)
{
  try {
    if (ch_ID == darkart::Channel::SUMCH_ID)
      return sumchannel;
    else {
      for (size_t i=0; i<channels.size(); i++) {
        if (ch_ID == channels[i].channel.channel_id())
          return channels[i];
      }
    }
    // couldn't find channel with desired channel ID
    throw std::invalid_argument("Requested channel could not be found! ");
  }
  catch (std::invalid_argument &err) {
    std::cout << err.what() << ch_ID <<std::endl;
  }
  throw std::runtime_error("");
}



#endif

