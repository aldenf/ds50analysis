/*

  Definition of the Channel class

  AFan 2013-07-01

*/


#ifndef darkart_Products_Channel_hh
#define darkart_Products_Channel_hh

#include <vector>
#include <string>
#include <iostream>

// Forward declarations
namespace darkart {
  struct Channel;
  typedef std::vector<Channel> ChannelVec;
}



struct darkart::Channel
{
  static const int INVALID_CHANNEL_ID = -1;
  static const int SUMCH_ID = -2;
  
  Channel():
    channelID(),
    board_id(-1),
    channel_num(-1),
    label(""),
    sample_bits(-1),
    sample_rate(0.),
    trigger_index(-1),
    nsamps(-1),
    saturated(false)
  { }

  struct ChannelID {
    int channel_id;
    int channel_index; 
  };

  ChannelID channelID;

  // All the identifier info we want to store for each channel
  int           board_id; //id of the board the channel is on
  int        channel_num; //number of the channel on the board
  std::string      label; //helpful label for each channel
  int      trigger_count; //used to check for ID_MISMATCH
  int        sample_bits; //digitizer resolution
  double     sample_rate; //samples per microsecond
  int      trigger_index; //sample at which trigger occurred
  int             nsamps; //number of samples in the waveform
  bool         saturated; //did the signal hit the max or min range of the digitizer?

  int channel_id() const { return channelID.channel_id; }
  double GetVerticalRange() const { return (1ULL<<sample_bits) - 1; }

  bool isValid() const {
    return (channel_id()!=-1 && sample_bits>-1 &&
            sample_rate>0 && trigger_index>-1 && nsamps>-1 && !saturated);
  }

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "Channel Info:"                                               << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << "Label:                      " << label                       << std::endl;
    std::cout << "Samples / microsec:         " << sample_rate                 << std::endl;
    std::cout << "Number of Samples:          " << nsamps                      << std::endl;
    std::cout << "Saturated:                  " << saturated                   << std::endl;

    if(verbosity > 0){
      std::cout << "Board ID:                   " << board_id                    << std::endl;
      std::cout << "Channel Number:             " << channel_num                 << std::endl;
      std::cout << "Trigger Count:              " << trigger_count               << std::endl;
      std::cout << "Resolution (bits):          " << sample_bits                 << std::endl;
      std::cout << "Trigger Index:              " << trigger_index               << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }
  
};


inline
bool operator==(darkart::Channel const& a, darkart::Channel const& b)
{
  return
    a.channelID.channel_id == b.channelID.channel_id &&
    a.channelID.channel_index == b.channelID.channel_index &&
    a.board_id      == b.board_id     &&
    a.channel_num   == b.channel_num;
}

#endif
