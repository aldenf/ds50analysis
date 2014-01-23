/*
  ROOT macro to play around with stuff
*/


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>
#include <tr1/cmath>

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TString.h"

#include "../darkart/Products/EventData.hh"
#include "analysis_lib.hh"

using namespace std;
using namespace darkart;

// Forward declaration
void LoopOverChain(TChain* tpc_chain, TString outFileName = "analysis.root");




void analysis() {

  TChain* tpc_chain = new TChain("treeBuilder/Events");
  string tpc_path = "/ds50/data/test_processing/darkart_release3/Run";
  std::vector<int> run_id_list;
  run_id_list.push_back(5370);
  run_id_list.push_back(5372);
  run_id_list.push_back(5373);
  /*
  run_id_list.push_back(5375);
  run_id_list.push_back(5376);
  run_id_list.push_back(5378);
  run_id_list.push_back(5379);
  run_id_list.push_back(5381);
  run_id_list.push_back(5382);
  run_id_list.push_back(5384);
  run_id_list.push_back(5385);
  run_id_list.push_back(5387);
  run_id_list.push_back(5388);
  run_id_list.push_back(5390);
  run_id_list.push_back(5391);
  run_id_list.push_back(5393);
  run_id_list.push_back(5396);
  run_id_list.push_back(5398);
  run_id_list.push_back(5399);
  */
 
  std::cout<<"WARNING: Database access disabled ! Make sure to check run list manually on the ELOG"<<std::endl;
  std::ostringstream os;
  for(vector<int>::const_iterator it = run_id_list.begin(); it != run_id_list.end(); it++)
    {
      os.str("");
      os << tpc_path << setw(6) << setfill('0') << *it
         << "/Run"<< setw(6) << setfill('0') << *it <<".root";
      std::cout<<"Adding file: "<<os.str().c_str()<<std::endl;
      tpc_chain->Add(os.str().c_str());
    }
    
  LoopOverChain(tpc_chain);
}


void analysis(TString Inputfilelist, TString outFileName) {

  TChain* tpc_chain = new TChain("treeBuilder/Events");

  std::cout<<"WARNING: Database access disabled ! Make sure to check run list manually on the ELOG"<<std::endl;

  Bool_t IsChained = AddFile2Chain(Inputfilelist, *tpc_chain);
  if (IsChained) LoopOverChain(tpc_chain, outFileName);
  else std::cout<<"Cannot chain files!! Please check input file list."<<std::endl;
}

void LoopOverChain(TChain* tpc_chain, TString outFileName)
{
  const Double_t t_drift_min = 10.; //[us]
  const Double_t t_drift_max = 376.; //[us]
  const Double_t t_drift_delta = 10; //[us]
  const Double_t electron_lifetime = 3330; //[us]
  const Int_t N_CHANNELS = 38;

  const Double_t s1_corr_limit = 2000.;
  const Int_t n_s1_slices = 20;
  const Int_t s1_slice_width = s1_corr_limit / n_s1_slices;

  const Double_t max_s2_top_frac_thresh = 0.4;

    
  Int_t tpc_events = tpc_chain->GetEntries();
  if (tpc_events == 0)
    {
      std::cout<<"No events loaded in TChain. Aborting..."<<std::endl;
      return;
    }
  std::cout << "Total number of events: "<<tpc_events << endl;
    
  EventData* event = NULL;
  tpc_chain->SetBranchAddress("EventData", &event);
    
  //CREATE OUTPUT FILES
  ofstream outfile;
  TString fOutliers(outFileName);
  fOutliers.Replace(fOutliers.Length()-5, 5, "_eventlist.txt"); //replace .root w/ _outliers.txt
  outfile.open(fOutliers.Data());
    
  //////////////////////////////////////////////////////////////////////////
  ///////////////////////     Declare histograms     ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  
  std::cout << "Saving output to "<<outFileName.Data()<<std::endl;
  TFile* f = new TFile(outFileName.Data(), "RECREATE");
  TH1F* t_drift_hist                = new TH1F("t_drift_hist", "Drift Time; t_drift [us]", 400, 0., 400.);
  TH1F* s1_startime_hist            = new TH1F("s1_startime_hist", "S1 Start Time; time [us]", 500, -1., 1.);
  TH1F* total_s1_hist               = new TH1F("total_s1_hist", "S1 Spectrum; S1 [p.e.]", 10000, 0, 10000);
  TH1F* total_s1_corr_hist          = new TH1F("total_s1_corr_hist", "S1 Spectrum (corrected for z-dependence); S1 [p.e.]",
                                               10000, 0, 10000);
  TH1F* total_s2_hist               = new TH1F("total_s2_hist", "S2 Spectrum; S2 [p.e.]", 5000, 0, 500000);
  TH1F* total_s2_corr_hist          = new TH1F("total_s2_corr_hist", "S2 Spectrum (corrected for z-dependence); S2 [p.e.]",
                                               5000, 0, 500000);
  TH1F* total_f90_hist              = new TH1F("total_f90_hist", "F90 Distribution; F90", 110, 0, 1.3);
  TH2F* total_s1_f90_hist           = new TH2F("total_s1_f90_hist", "F90 vs S1; S1 [p.e.]; F90", 10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_hist      = new TH2F("total_s1_corr_f90_hist", "F90 vs S1 (corrected for z-dependence); S1 [p.e.]; F90",
                                               10000, 0, 10000, 130, 0, 1.3);
  TH2F* s2_s1_corr_hist             = new TH2F("s2_s1_corr_hist", "S2 vs S1 (both corr. for z-dep.); S1 [p.e.]; S2/S1",
                                               1000, 0, 6000, 1000, 0, 200);
  TH2F* s2_s1_corr_xsatcut_hist     = new TH2F("s2_s1_corr_xsatcut_hist", "S2 vs S1 (both corr. for z-dep.), no sat. cut; S1 [p.e.]; S2/S1",
                                               1000, 0, 6000, 1000, 0, 200);
  TH2F* s2_s1_corr_saturated_hist   = new TH2F("s2_s1_corr_saturated_hist", "S2 vs S1 (both corr. for z-dep.), ADC sat. only; S1 [p.e.]; S2/S1",
                                               1000, 0, 6000, 1000, 0, 200);


  TH2F* max_s2_top_s1_corr_hist     = new TH2F("max_s2_top_s1_corr_hist", "max S2 / top S2 vs S1_corr; S1 [p.e.]; max S2 / top S2",
                                               400, 0., 5000, 200, 0, 1.);
  TH2F* max_s2_bot_s1_corr_hist     = new TH2F("max_s2_bot_s1_corr_hist", "max S2 / bot S2 vs S1 (corr.); S1 [p.e.]; max S2 / bot S2",
                                               400, 0., 5000, 200, 0, 2.);
  
  TH2F* xy_hist                    = new TH2F("xy_hist", "x-y; x [arb]; y [arb]", 500, -4, 4, 500, -4, 4);
  
  TH2F* s2_length_s1_corr_hist     = new TH2F("s2_length_s1_corr_hist", "S2 pulse length vs S1 (corr.); S1 [p.e.]; Length [us]",
                                              400, 0., 5000, 200, 0, 400);
  TH2F* s2fullfixed_s1_corr_hist   = new TH2F("s2fullfixed_s1_corr_hist", "S2 (full) / S2 (fixed) vs S1 (z-corr.); S1 [p.e.]; S2 (full) / S1 (fixed)",
                                              400, 0, 5000, 1000, 0, 2);


  
  TH1F** max_s2_top_hists = new TH1F*[n_s1_slices];
  TH2F** xy_hists = new TH2F*[n_s1_slices];
  TH2F** xy_under_pmt_hists = new TH2F*[n_s1_slices];
  TH2F** xy_xunder_pmt_hists = new TH2F*[n_s1_slices];
  TH1F** total_s2_corr_hists = new TH1F*[n_s1_slices];
  TH1F** total_s2_corr_under_pmt_hists = new TH1F*[n_s1_slices];
  TH1F** total_s2_corr_xunder_pmt_hists = new TH1F*[n_s1_slices];
  TH1F** bot_s2_under_pmt_hists = new TH1F*[n_s1_slices];
  TH1F** bot_s2_xunder_pmt_hists = new TH1F*[n_s1_slices];
  TH1F** bot_s2_corr_hists = new TH1F*[n_s1_slices];
  TH1F** s2_length_s1_corr_hists = new TH1F*[n_s1_slices];
  TH1F** s2fullfixed_s1_corr_hists = new TH1F*[n_s1_slices];

  TH1F*** total_s2_corr_hists_chans = new TH1F**[n_s1_slices];
  TH1F*** max_s2_top_hists_chans = new TH1F**[n_s1_slices];
  TH1F*** max_s2_bot_hists_chans = new TH1F**[n_s1_slices];
  TH1F*** total_s2_corr_hists_radii = new TH1F**[n_s1_slices];
  TH1F*** bot_s2_corr_hists_radii = new TH1F**[n_s1_slices];
  
  
  
  
  for (int i=0; i<n_s1_slices; i++) {

    std::ostringstream name, title;
    name<<"max_s2_top_hist_"<<i;
    title<<"max S2 / top S2, "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; max_s2/top_s2";
    max_s2_top_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 500, 0, 1.);

    name.str(""); title.str("");
    name<<"xy_hist_"<<i;
    title<<"x-y, "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; x [arb]; y [arb";
    xy_hists[i] = new TH2F(name.str().c_str(), title.str().c_str(), 500, -4, 4, 500, -4, 4);

    name.str(""); title.str("");
    name<<"xy_under_pmt_hist_"<<i;
    title<<"x-y, max_s2_top_frac > "<<max_s2_top_frac_thresh<<", "
         <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; x [arb]; y [arb";
    xy_under_pmt_hists[i] = new TH2F(name.str().c_str(), title.str().c_str(), 500, -4, 4, 500, -4, 4);

    name.str(""); title.str("");
    name<<"xy_xunder_pmt_hist_"<<i;
    title<<"x-y, max_s2_top_frac < "<<max_s2_top_frac_thresh<<", "
         <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; x [arb]; y [arb]";
    xy_xunder_pmt_hists[i] = new TH2F(name.str().c_str(), title.str().c_str(), 500, -4, 4, 500, -4, 4);

    name.str(""); title.str("");
    name<<"total_s2_corr_hist_"<<i;
    title<<"S2 (corr.), "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
    total_s2_corr_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 200000);
    
    name.str(""); title.str("");
    name<<"total_s2_corr_under_pmt_hist_"<<i;
    title<<"S2 (corr.), max_s2_top_frac > "<<max_s2_top_frac_thresh<<", "
         <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
    total_s2_corr_under_pmt_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 200000);

    name.str(""); title.str("");
    name<<"total_s2_corr_xunder_pmt_hist_"<<i;
    title<<"S2 (corr.), max_s2_top_frac < "<<max_s2_top_frac_thresh<<", "
         <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
    total_s2_corr_xunder_pmt_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 200000);

    name.str(""); title.str("");
    name<<"bot_s2_under_pmt_hist_"<<i;
    title<<"bot S2 (corr.), max_s2_top_frac > "<<max_s2_top_frac_thresh<<", "
         <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
    bot_s2_under_pmt_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 100000);

    name.str(""); title.str("");
    name<<"bot_s2_xunder_pmt_hist_"<<i;
    title<<"bot S2 (corr.), max_s2_top_frac < "<<max_s2_top_frac_thresh<<", "
         <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
    bot_s2_xunder_pmt_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 100000);

    name.str(""); title.str("");
    name<<"bot_s2_corr_hist_"<<i;
    title<<"Bottom S2 (corr.), "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
    bot_s2_corr_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 100000);


    name.str(""); title.str("");
    name<<"s2_length_s1_corr_hist_"<<i;
    title<<"S2 pulse length, "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; length [us]";
    s2_length_s1_corr_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 400);

    name.str(""); title.str("");
    name<<"s2fullfixed_s1_corr_hist_"<<i;
    title<<"S2 (full) / S2 (fixed), "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 (full) / S2 (fixed)";
    s2fullfixed_s1_corr_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 2);

    
    total_s2_corr_hists_chans[i] = new TH1F*[N_CHANNELS];
    max_s2_top_hists_chans[i] = new TH1F*[N_CHANNELS];
    max_s2_bot_hists_chans[i] = new TH1F*[N_CHANNELS];
    for (int j=0; j<N_CHANNELS; j++) {
      name.str(""); title.str("");
      name<<"total_s2_corr_under_pmt_hist_"<<i<<"_ch"<<j;
      title<<"S2 (corr.), max_s2_top_frac > "<<max_s2_top_frac_thresh<<", max S2 ch"<<j<<", "
           <<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
      total_s2_corr_hists_chans[i][j] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 200000);

      name.str(""); title.str("");
      name<<"max_s2_top_hist_"<<i<<"_ch"<<j;
      title<<"max S2 / top S2, ch"<<j<<", "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; max_s2/top_s2";
      max_s2_top_hists_chans[i][j] = new TH1F(name.str().c_str(), title.str().c_str(), 500, 0, 1.);

      name.str(""); title.str("");
      name<<"max_s2_bot_hist_"<<i<<"_ch"<<j;
      title<<"max S2 / bot S2, ch"<<j<<", "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; max_s2/bot_s2";
      max_s2_bot_hists_chans[i][j] = new TH1F(name.str().c_str(), title.str().c_str(), 500, 0, 2.);
    }

    total_s2_corr_hists_radii[i] = new TH1F*[4];
    bot_s2_corr_hists_radii[i] = new TH1F*[4];
    for (int j=0; j<4; j++) {
      name.str(""); title.str("");
      name<<"total_s2_corr_hist_"<<i<<"_radius_"<<j;
      title<<"S2 (corr.), radius "<<j<<", "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
      total_s2_corr_hists_radii[i][j] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 200000);

      name.str(""); title.str("");
      name<<"bot_s2_corr_hist_"<<i<<"_radius_"<<j;
      title<<"Bottom S2 (corr.), radius "<<j<<", "<<i*s1_slice_width<<" < S1 corr. < "<<(i+1)*s1_slice_width<<"; S2 [pe]";
      bot_s2_corr_hists_radii[i][j] = new TH1F(name.str().c_str(), title.str().c_str(), 1000, 0, 100000);
    }
  }


  //////////////////////////////////////////////////////////////////////////
  /////////////////     BEGIN LOOP OVER EVENTS       ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  //tpc_events = tpc_events - 50; // Skip last few events because end of some (very few) runs are problematic.
  for (int n = 0; n < tpc_events; n++)
    //for (int n=0; n<10000; n++)
    {
        
      //Load the event
      tpc_chain->GetEntry(n);
        
      if ( n % 10000 == 0)
        std::cout<<"Processing Event: "<<n<<"/"<<tpc_events<<std::endl;

      /////////////////////////
      //  APPLY BASIC CUTS   //
      /////////////////////////

      // Check for expected number of channels
      if ((int)event->channels.size() != N_CHANNELS){
        cout << "Event=" << event->event_info.event_id<<" has LOWER # of Channels"<<endl;
        continue;
      }
      
      //Make sure the event is not saturated
      //if (event->event_info.saturated == true)
      //  continue;
        
      //Make sure the baseline was found on the sum channel
      if (event->sumchannel.baseline.found_baseline == false)
        continue;
        
      //PULSE IDENTIFICATION
      Int_t n_phys_pulses = -1, s1_pulse_id = -1, s2_pulse_id = -1;
      ds50analysis::identify_pulses(event, n_phys_pulses, s1_pulse_id, s2_pulse_id, t_drift_max, t_drift_delta);

      
      //Make sure there are 2 pulses
      if (n_phys_pulses != 2)
        continue;
        
      //CALCULATE PARAMETERS
      Double_t total_s1 = event->pulses[s1_pulse_id].param.fixed_int1;
        
      Double_t s1_start_time = event->sumchannel.pulses[s1_pulse_id].pulse.start_time;
      Double_t s2_start_time = event->sumchannel.pulses[s2_pulse_id].pulse.start_time;
      Double_t t_drift = s2_start_time - s1_start_time;

      Double_t s2_pulse_length = event->sumchannel.pulses[s2_pulse_id].pulse.end_time - s2_start_time;
        
      Double_t total_s1_corr = total_s1*ds50analysis::s1_corr_factor(t_drift_max, t_drift);
      Double_t total_s1_90 = event->pulses[s1_pulse_id].param.f90*event->pulses[s1_pulse_id].param.npe;
      Double_t total_f90 = total_s1_90/total_s1;
        
      Double_t total_s2 = event->pulses[s2_pulse_id].param.fixed_int2;
      Double_t total_s2_corr = total_s2/TMath::Exp(-t_drift/electron_lifetime);
      Double_t total_s2_full = event->pulses[s2_pulse_id].param.integral;
      Double_t s2full_over_s2fixed = total_s2_full / total_s2;
      
      Double_t s2_over_s1_corr = total_s2_corr / total_s1_corr;

      Int_t max_s1_chan = -1, max_s2_chan = -1;
      Double_t max_s1 = -1., max_s2 = -1.;
      ds50analysis::max_s1_s2(event, s1_pulse_id, s2_pulse_id, max_s1_chan, max_s1, max_s2_chan, max_s2);

      Double_t top_s1=0, bot_s1=0, top_s2=0, bot_s2=0;
      ds50analysis::top_bot_s1_s2(event, s1_pulse_id, s2_pulse_id, top_s1, bot_s1, top_s2, bot_s2);

      Double_t bot_s2_corr = bot_s2 / TMath::Exp(-t_drift/electron_lifetime);
      
      //Double_t max_s1_frac = max_s1/total_s1;
      //Double_t max_s2_frac = max_s2/total_s2;
      Double_t max_s2_top_frac = max_s2 / top_s2;
      Double_t max_s2_bot_frac = max_s2 / bot_s2;

      Double_t max_s2_chan_radius = TMath::Sqrt(event->getChannelByID(max_s2_chan).pmt.photocathode_x*
                                                event->getChannelByID(max_s2_chan).pmt.photocathode_x+
                                                event->getChannelByID(max_s2_chan).pmt.photocathode_y*
                                                event->getChannelByID(max_s2_chan).pmt.photocathode_y);
      
      Double_t xpos = event->pulses[s2_pulse_id].position.bary_x;
      Double_t ypos = event->pulses[s2_pulse_id].position.bary_y;

      
      s1_startime_hist->Fill(s1_start_time);
      
      ///////////////////////////
      //  APPLY ANALYSIS CUTS  //
      ///////////////////////////
      
      //Make sure the S1 pulse is where you expect in the trigger window
      if (s1_start_time < -0.25 || s1_start_time > -0.18 )
        continue;

      t_drift_hist->Fill(t_drift);
      
      //Remove very slow events (triggered on S2)
      if (total_f90 < 0.1)
        continue;
        
      //Remove events concentrated on a single PMT
      if (max_s1/total_s1 > 0.4)
        continue;
        
      //Remove events near grid or cathode
      if (t_drift < t_drift_min || t_drift > t_drift_max - t_drift_delta)
        continue;
      
      // Cut on large max_s1_frac. Threshold is defined bin by bin.
      if (ds50analysis::large_max_s1_frac(total_s1, t_drift, max_s1))
        continue;


      s2_s1_corr_xsatcut_hist->Fill(total_s1_corr, s2_over_s1_corr);
      
      //Make sure the event is not saturated
      if (event->event_info.saturated == true) {
        s2_s1_corr_saturated_hist->Fill(total_s1_corr, s2_over_s1_corr);
        continue;
      }




     

      

      ////////////////////////
      // Fill eventlist file.
      // darkart::InEventlist module expects first 2 entries of each column to be run_id & event_id.
      // More entries can be included if other modules (e.g. AverageWaveforms) needs it.
      if ( (t_drift > 30. && t_drift < 40.) &&
           (total_s1_corr > 0 && total_s1_corr < 2000) /*&& (max_s2_chan == 30)*/ )
        {
          outfile << event->event_info.run_id <<" "<< event->event_info.event_id <<" "<< s2_start_time<<" "
                  << total_s1_corr<<" "<<max_s2_chan
                  << std::endl;
        }
      ////////////////////////

      
      ///////////////////////
      //  FILL HISTOGRAMS  //
      ///////////////////////
      
      
      total_s1_hist               ->Fill(total_s1);
      total_s1_corr_hist          ->Fill(total_s1_corr);
      total_s2_hist               ->Fill(total_s2);
      total_s2_corr_hist          ->Fill(total_s2_corr);
      total_f90_hist              ->Fill(total_f90);
      total_s1_f90_hist           ->Fill(total_s1, total_f90);
      total_s1_corr_f90_hist      ->Fill(total_s1_corr, total_f90);
      s2_s1_corr_hist             ->Fill(total_s1_corr, s2_over_s1_corr);

      max_s2_top_s1_corr_hist->Fill(total_s1_corr, max_s2_top_frac);
      max_s2_bot_s1_corr_hist->Fill(total_s1_corr, max_s2_bot_frac);
      
      xy_hist->Fill(xpos, ypos);
      
      s2_length_s1_corr_hist->Fill(total_s1_corr, s2_pulse_length);
      s2fullfixed_s1_corr_hist->Fill(total_s1_corr, s2full_over_s2fixed);



      for (int i=0; i<n_s1_slices; i++) {
        if (total_s1_corr > i*s1_slice_width && total_s1_corr <= (i+1)*s1_slice_width) {
          max_s2_top_hists[i]->Fill(max_s2_top_frac);
          xy_hists[i]->Fill(xpos, ypos);
          total_s2_corr_hists[i]->Fill(total_s2_corr);
          max_s2_top_hists_chans[i][max_s2_chan]->Fill(max_s2_top_frac);
          max_s2_bot_hists_chans[i][max_s2_chan]->Fill(max_s2_bot_frac);
          bot_s2_corr_hists[i]->Fill(bot_s2_corr);
          s2_length_s1_corr_hists[i]->Fill(s2_pulse_length);
          s2fullfixed_s1_corr_hists[i]->Fill(s2full_over_s2fixed);
          if (max_s2_top_frac > max_s2_top_frac_thresh) {
            xy_under_pmt_hists[i]->Fill(xpos, ypos);
            total_s2_corr_under_pmt_hists[i]->Fill(total_s2_corr);
            bot_s2_under_pmt_hists[i]->Fill(bot_s2);
            total_s2_corr_hists_chans[i][max_s2_chan]->Fill(total_s2_corr);
          }
          else {
            xy_xunder_pmt_hists[i]->Fill(xpos, ypos);
            total_s2_corr_xunder_pmt_hists[i]->Fill(total_s2_corr);
            bot_s2_xunder_pmt_hists[i]->Fill(bot_s2);
          }
          if (max_s2_chan_radius >= 0 && max_s2_chan_radius < 1.) {
            total_s2_corr_hists_radii[i][0]->Fill(total_s2_corr);
            bot_s2_corr_hists_radii[i][0]->Fill(bot_s2_corr);
          }
          else if (max_s2_chan_radius > 1. && max_s2_chan_radius < 3.) {
            total_s2_corr_hists_radii[i][1]->Fill(total_s2_corr);
            bot_s2_corr_hists_radii[i][1]->Fill(bot_s2_corr);
          }
          else if (max_s2_chan_radius > 3. && max_s2_chan_radius < 3.5) {
            total_s2_corr_hists_radii[i][2]->Fill(total_s2_corr);
            bot_s2_corr_hists_radii[i][2]->Fill(bot_s2_corr);
          }
          else {
            total_s2_corr_hists_radii[i][3]->Fill(total_s2_corr);
            bot_s2_corr_hists_radii[i][3]->Fill(bot_s2_corr);
          }
        }

      }





      
    }//End loop over events

  total_s1_f90_hist->SetDrawOption("COLZ");
  total_s1_corr_f90_hist->SetDrawOption("COLZ");

  /////////////////////////////////////////////////////////////////////////
  /////////////////     WRITE HISTOGRAMS TO FILE       ////////////////////
  /////////////////////////////////////////////////////////////////////////
  t_drift_hist->Write();
  s1_startime_hist->Write();
  total_s1_hist->Write();
  total_s1_corr_hist->Write();
  total_s2_hist->Write();
  total_s2_corr_hist->Write();
  total_f90_hist->Write();
  total_s1_f90_hist->Write();
  total_s1_corr_f90_hist->Write();
  s2_s1_corr_hist->Write();
  s2_s1_corr_xsatcut_hist->Write();
  s2_s1_corr_saturated_hist->Write();

  max_s2_top_s1_corr_hist->Write();
  max_s2_bot_s1_corr_hist->Write();
  xy_hist->Write();
  s2_length_s1_corr_hist->Write();
  s2fullfixed_s1_corr_hist->Write();

  for (int i=0; i<n_s1_slices; i++) {

    std::ostringstream dirname;
    dirname<<"s1_corr_slice_"<<i;
    f->mkdir(dirname.str().c_str());
    f->cd(dirname.str().c_str());
    
    
    max_s2_top_hists[i]->Write();
    xy_hists[i]->Write();
    xy_under_pmt_hists[i]->Write();
    xy_xunder_pmt_hists[i]->Write();
    total_s2_corr_hists[i]->Write();
    total_s2_corr_under_pmt_hists[i]->Write();
    total_s2_corr_xunder_pmt_hists[i]->Write();
    bot_s2_under_pmt_hists[i]->Write();
    bot_s2_xunder_pmt_hists[i]->Write();
    bot_s2_corr_hists[i]->Write();
    s2_length_s1_corr_hists[i]->Write();
    s2fullfixed_s1_corr_hists[i]->Write();

    dirname.str("");
    dirname<<"s1_corr_slice_"<<i<<"/channels";
    f->mkdir(dirname.str().c_str());
    f->cd(dirname.str().c_str());
    for (int j=0; j<N_CHANNELS; j++) {
      total_s2_corr_hists_chans[i][j]->Write();
      max_s2_top_hists_chans[i][j]->Write();
      max_s2_bot_hists_chans[i][j]->Write();
    }

    dirname.str("");
    dirname<<"s1_corr_slice_"<<i<<"/radii";
    f->mkdir(dirname.str().c_str());
    f->cd(dirname.str().c_str());
    for (int j=0; j<4; j++) {
      total_s2_corr_hists_radii[i][j]->Write();
      bot_s2_corr_hists_radii[i][j]->Write();
    }
  }
  
  //new TCanvas();
  //total_s1_hist->Draw("");
  //new TCanvas();
  //total_s1_corr_hist->Draw("");
  f->cd();
  f->Close();
  outfile.close();


}




#ifndef __CINT__
int main(int argc, char **argv) {
  if ( argc == 1 ) {
    std::cout << "\n==========> analysis <=============" << std::endl;
    analysis();
  } else if ( argc == 3 ) {
    std::cout << "\n==========> analysis with file list & outputname <=============" << std::endl;
    analysis(argv[1], argv[2]);
  } else {
    std::cout << "Usage:" <<argc<< std::endl;
    std::cout << "./analysis" << std::endl;
    std::cout << "OR ./analysis filelist.txt output.root" << std::endl;
    return 0;
  }


  std::cout << "==> Application finished." << std::endl;

  return 0;
}
#endif /* __CINT __ */

