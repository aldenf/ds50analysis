/*

  Suite of functions to be used in analysis macros.

*/

#ifndef analysis_g2proposal_lib_hh
#define analysis_g2proposal_lib_hh

#include <vector>
#include <string>
#include <fstream>


#include "TFile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"

#include "../darkart/Products/EventData.hh"


using namespace std;
using namespace darkart;


namespace ds50analysis {
  
  // Geometrical correction factor for S1.
  Double_t s1_corr_factor(Double_t t_drift_max, Double_t t_drift);

  // Identify the number of physical pulses and, if possible, S1 pulse ID and S2 pulse ID.
  void identify_pulses(EventData* event,
                       Int_t & n_phys_pulses, Int_t & s1_pulse_id, Int_t & s2_pulse_id,
                       Double_t t_drift_max, Double_t t_drift_delta);


  // max_chan calculation in DarkArt uses npe and not npe_fixed
  // --> must recalculate here using fixed_int params
  void max_s1_s2(EventData* event, Int_t const s1_pulse_id, Int_t const s2_pulse_id,
                 Int_t & max_s1_chan, Double_t & max_s1,
                 Int_t & max_s2_chan, Double_t & max_s2);

  
  // xy barycenter is not yet done in DarkArt, so do it here.
  void barycenter(EventData* event, Int_t const s2_pulse_id,
                  Double_t & xpos, Double_t & ypos,
                  Double_t const* pmt_xpos, Double_t const* pmt_ypos,
                  Int_t const* topChannels);

  // Cut on large max_s1_frac. Threshold is defined bin by bin.
  Bool_t large_max_s1_frac(Double_t const total_s1, Double_t const t_drift, Double_t const max_s1);


  

  
}








Double_t ds50analysis::s1_corr_factor(Double_t t_drift_max, Double_t t_drift)
{
  Double_t t_drift_ratio = t_drift/(0.5*t_drift_max); // note normalization is to 0.5*t_drift_max
  // looked at Kr peaks in 25us drift-time windows (Run5330+5340), and fit these to [0]*t_drift_ratio^2 + [1]*t_drift_ratio + [2].
  Double_t fit_par0 = 4.458;
  Double_t fit_par1 = 9.805;
  Double_t fit_par2 = 281.2;
  // normalizing all points on fitted curve to expected Kr peak at t_drift_max/2
  Double_t exp_Kr_peak_at_half_t_drift_max = fit_par0 + fit_par1 + fit_par2;
  Double_t exp_Kr_peak_at_t_drift = fit_par0*t_drift_ratio*t_drift_ratio + fit_par1*t_drift_ratio + fit_par2;
  return exp_Kr_peak_at_half_t_drift_max/exp_Kr_peak_at_t_drift; // s1 correction factor
}



void ds50analysis::identify_pulses(EventData* event,
                                   Int_t & n_phys_pulses, Int_t & s1_pulse_id, Int_t & s2_pulse_id,
                                   Double_t t_drift_max, Double_t t_drift_delta)
{
  if (event->sumchannel.pulses.size() == 0)
    {
      n_phys_pulses = 0;
    }
  else if (event->sumchannel.pulses.size() == 1)
    {
      //Assume pulse is S1 ... for now
      n_phys_pulses = 1;
      s1_pulse_id = 0;
    }
  else if (event->sumchannel.pulses.size() == 2)
    {
      //Assume first pulse is S1 and second is S2 ... for now
      n_phys_pulses = 2;
      s1_pulse_id = 0;
      s2_pulse_id = 1;
    }
  else if (event->sumchannel.pulses.size() == 3
           && std::fabs(event->sumchannel.pulses[2].pulse.start_time
                        - event->sumchannel.pulses[1].pulse.start_time
                        - t_drift_max) < t_drift_delta)
    {
      //Assume first pulse is S1, second is S2 and third is S3 ... for now
      n_phys_pulses = 2;
      s1_pulse_id = 0;
      s2_pulse_id  = 1;
    }
  else
    { //We don't know how many physical pulses - just set to total number of pulses for now
      n_phys_pulses = event->sumchannel.pulses.size();
    }
}

void ds50analysis::max_s1_s2(EventData* event, Int_t const s1_pulse_id, Int_t const s2_pulse_id,
                             Int_t & max_s1_chan, Double_t & max_s1,
                             Int_t & max_s2_chan, Double_t & max_s2)
{
  const Int_t nchans = event->channels.size();
  for (Int_t ch = 0; ch < nchans; ch++)
    {
      ChannelData const& channel = event->getChannelByID(ch);
      Double_t s1 = -channel.pulses[s1_pulse_id].param.fixed_int1 / channel.pmt.spe_mean;
      Double_t s2 = -channel.pulses[s2_pulse_id].param.fixed_int2 / channel.pmt.spe_mean;
      if (s1 > max_s1)
        {
          max_s1 = s1;
          max_s1_chan = ch;
        }
      if (s2 > max_s2)
        {
          max_s2 = s2;
          max_s2_chan = ch;
        }
    }
  
}




void ds50analysis::barycenter(EventData* event, Int_t const s2_pulse_id,
                              Double_t & xpos, Double_t & ypos,
                              Double_t const* pmt_xpos, Double_t const* pmt_ypos,
                              Int_t const* topChannels)
{
  xpos = 0;
  ypos = 0;
  Double_t total_s2_top = 0;
  const Int_t nchans = event->channels.size();
  for (Int_t ch = 0; ch < nchans; ch++)
    {
      ChannelData const& channel = event->getChannelByID(ch);
      if (channel.channel.channel_id()>=19)
        {
          total_s2_top += -channel.pulses[s2_pulse_id].param.fixed_int2 / channel.pmt.spe_mean;
          for (Int_t j=0; j<19; j++)
            {
              if (topChannels[j] == channel.channel.channel_id())
                {
                  xpos += pmt_xpos[j] * -channel.pulses[s2_pulse_id].param.fixed_int2 / channel.pmt.spe_mean;
                  ypos += pmt_ypos[j] * -channel.pulses[s2_pulse_id].param.fixed_int2 / channel.pmt.spe_mean;
                    
                  break;
                }
            }
        }

    }//end loop over channels
  xpos /= total_s2_top;
  ypos /= total_s2_top;

}



Bool_t ds50analysis::large_max_s1_frac(Double_t const total_s1, Double_t const t_drift, Double_t const max_s1)
{
  Bool_t result = false;
  
  // Fill the histograms
  if (total_s1>=40. && total_s1<50. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=50. && total_s1<60. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.18)
    result = true;
  else if (total_s1>=60. && total_s1<70. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.16)
    result = true;
  else if (total_s1>=70. && total_s1<80. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.16)
    result = true;
  else if (total_s1>=80. && total_s1<90. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.15)
    result = true;
  else if (total_s1>=90. && total_s1<100. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.15)
    result = true;
  else if (total_s1>=100. && total_s1<110. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.15)
    result = true;
  else if (total_s1>=110. && total_s1<120. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.15)
    result = true;
  else if (total_s1>=120. && t_drift>=50. && t_drift<=300. && max_s1/total_s1 > 0.15)
    result = true;
        
  else if (total_s1>=40. && total_s1<50. && t_drift<50. && max_s1/total_s1 > 0.25)
    result = true;
  else if (total_s1>=50. && total_s1<60. && t_drift<50. && max_s1/total_s1 > 0.25)
    result = true;
  else if (total_s1>=60. && total_s1<70. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=70. && total_s1<80. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=80. && total_s1<90. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=90. && total_s1<100. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=100. && total_s1<110. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=110. && total_s1<120. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
  else if (total_s1>=120. && t_drift<50. && max_s1/total_s1 > 0.2)
    result = true;
        
  else if (total_s1>=40. && total_s1<50. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=50. && total_s1<60. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=60. && total_s1<70. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=70. && total_s1<80. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=80. && total_s1<90. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=90. && total_s1<100. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=100. && total_s1<110. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=110. && total_s1<120. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;
  else if (total_s1>=120. && t_drift>300. && max_s1/total_s1 > 0.4)
    result = true;

  return result;

}

  
Bool_t AddFile2Chain(TString Inputfilelist, TChain &chain){
  Bool_t IsChained(false);
  ifstream inputStream(Inputfilelist.Data());
  if (!inputStream.is_open()) {
    cout << "can not open list file"<< endl;
    return false;
  }
  cout<<"Open file list: "<<Inputfilelist.Data()<<endl;

  char line[512];
  for (; inputStream.good();) {
    inputStream.getline(line, 512);
    if (!inputStream.good()) continue;

    TFile *ftmp = new TFile(line);
    //----------
    if (!(ftmp->IsOpen())){
      cout << line << " open failed ! not chained" << endl;
      continue;
    }
    if (ftmp->IsZombie()) {
      cout << "sth. very wrong with " << line << ", not chained " << endl;
      continue;
    }
    if (ftmp->TestBit(1024)) {
      cout << "revocer procedure applied to " << line << endl;
      continue;
    }
    //--------------------------
    if (ftmp && ftmp->IsOpen() && ftmp->GetNkeys()) {
      cout << "add file " << line << endl;
      chain.Add(line);
      IsChained=true;
    } else {
      cout << " cannot open file " << line << endl;
    }
    delete ftmp;
  }
  return IsChained;
}
#endif
