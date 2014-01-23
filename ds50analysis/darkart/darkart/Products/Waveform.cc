#include "darkart/Products/Waveform.hh"
#include "art/Utilities/Exception.h"
#include <limits>


bool darkart::Waveform::isValid() const
{
  return (!wave.empty() && trigger_index > -1 && sample_rate > 0);
}

double darkart::Waveform::SampleToTime(int sample) const
{
  if (!isValid()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return ((double)(sample - trigger_index)) / sample_rate;
}

int darkart::Waveform::TimeToSample(double time, bool checkrange) const
{
  if (!isValid()) {
    return -1;
  }
  
  int nsamps = wave.size();
  int samp = (int)(time*sample_rate + trigger_index+0.001); //add 0.001 to avoid rounding errors
  if (checkrange) {
    if (samp < 0) samp = 0;
    if (samp > nsamps -1) samp = nsamps - 1;
  }
  return samp;
}

