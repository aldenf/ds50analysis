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

#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TNtuple.h"

#include "../darkart/Products/EventData.hh"
#include "analysis_jan2014note_lib.hh"

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
  fOutliers.Replace(fOutliers.Length()-5, 5, "_outliers.txt"); //replace .root w/ _outliers.txt
  outfile.open(fOutliers.Data());
    
  //////////////////////////////////////////////////////////////////////////
  ///////////////////////     Declare histograms     ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  
  std::cout << "Saving output to "<<outFileName.Data()<<std::endl;
  TFile* f = new TFile(outFileName.Data(), "RECREATE");
  TH1F* t_drift_hist                = new TH1F("t_drift_hist", "Drift Time", 400, 0., 400.);
  TH1F* s1_startime_hist            = new TH1F("s1_startime_hist", "Drift Time", 500, -1., 1.);
  TH1F* total_s1_hist               = new TH1F("total_s1_hist", "S1 Spectrum", 10000, 0, 10000);
  TH1F* total_s1_corr_hist          = new TH1F("total_s1_corr_hist", "S1 Spectrum (corrected for z-dependence)", 10000, 0, 10000);
  TH1F* total_s2_hist               = new TH1F("total_s2_hist", "S2 Spectrum", 5000, 0, 500000);
  TH1F* total_s2_corr_hist          = new TH1F("total_s2_corr_hist", "S2 Spectrum (corrected for z-dependence)", 5000, 0, 500000);
  TH1F* total_f90_hist              = new TH1F("total_f90_hist", "F90 Distribution", 110, 0, 1.3);
  TH2F* total_s1_f90_hist           = new TH2F("total_s1_f90_hist", "F90 vs S1; S1 [p.e.]; F90", 10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_hist      = new TH2F("total_s1_corr_f90_hist", "F90 vs S1 (corrected for z-dependence); S1 [p.e.]; F90",
                                               10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f20_hist      = new TH2F("total_s1_corr_f20_hist", "F20 vs S1 (corrected for z-dependence); S1 [p.e.]; F90",
                                               10000, 0, 10000, 130, 0, 1.3);
  TH2F* max_frac_s1_hist            = new TH2F("max_frac_s1_hist", "Maximum S1 fraction vs S1",
                                               10000, 0, 10000, 100, 0, 1.0);
  TH2F* max_frac_max_peak_time_hist = new TH2F("max_frac_max_peak_time_hist", "Max Peak Arrival Time vs Maximum S1 fraction vs ",
                                               10000, 0, 1.0, 800, 0, 8.0);
  TH2F* s1_max_peak_time_hist       = new TH2F("s1_max_peak_time_hist", "Max Peak Arrival Time vs Maximum S1",
                                               10000, 0, 10000, 100, 0, 8.0);
  TH2F* t_drift_s1_corr_hist        = new TH2F("t_drift_s1_corr_hist", "Drift time vs S1 (corrected for z-dependence)",
                                               1000, 0, 6000, 400, 0., 400.);
  TH2F* s2_s1_corr_hist             = new TH2F("s2_s1_corr_hist", "S2 vs S1 (corrected for z-dependence)",
                                               1000, 0, 6000, 5000, 0, 50000);
  TH2F* logs2s1_corr_s1_corr_hist   = new TH2F("s2s1_corr_s1_corr_hist", "Log(S2/S1) vs S1 S1 (corrected for z-dependence)",
                                               1000, 0, 6000, 200, -1, 3);
  TH1F* total_livetime              = new TH1F("total_livetime", "total_livetime", 1, 0, 1);
  TH1F* max_s2_chan_hist            = new TH1F("max_s2_chan_hist", "Max S2 channels", 100, 0, 50);

  const int n_hists = 40;
  TH2F** logs2s1_corr_f90_hist = new TH2F*[n_hists];
  TH2F** max_frac_t_drift_hist = new TH2F*[n_hists];
  for (int i = 0; i < n_hists; i++)
    {
      std::ostringstream name, title;
      name<<"logs2s1_f90_hist_"<<i;
      title<<"Log(S2/S1) vs F90 for S1(corr) above "<<i*5<<" p.e.; F90; Log(S2/S1)";
      logs2s1_corr_f90_hist[i] = new TH2F(name.str().c_str(), title.str().c_str(), 120, 0.0, 1.2, 200, -1, 3);

      name.str("");title.str("");
      name<<"max_frac_t_drift_hist_"<<i;
      title<<"Max S1 Fraction vs Drift time for "<<(i)*5<<" < S1 [p.e.] < "<<(i+1)*5;
      max_frac_t_drift_hist[i] = new TH2F(name.str().c_str(), title.str().c_str(), 400, 0, 400, 100, 0, 1.0);
    }
    
  Double_t LivetimeTotal(0.), InhibitTimeTotal(0.);


  //////////////////////////////////////////////////////////////////////////
  /////////////////     BEGIN LOOP OVER EVENTS       ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  tpc_events = tpc_events - 50; // Skip last few events because end of some (very few) runs are problematic.
  for (int n = 0; n < tpc_events; n++)
    {
        
      //Load the event
      tpc_chain->GetEntry(n);
        
      if ( n % 10000 == 0)
        std::cout<<"Processing Event: "<<n<<"/"<<tpc_events<<std::endl;
        
      if(event->event_info.live_time_20ns*20.*1.e-9 < 1.)
        { //long Livetime event cut
          LivetimeTotal+=event->event_info.live_time_20ns*20.*1.e-9; // in second
          InhibitTimeTotal+=event->event_info.incremental_inhibit_time_20ns*20.*1.e-9; // in second
        }// This should be before any event cuts!!
        


      /////////////////////////
      //  APPLY BASIC CUTS   //
      /////////////////////////

      // Check for expected number of channels
      if ((int)event->channels.size() != N_CHANNELS){
        cout << "Event=" << event->event_info.event_id<<" has LOWER # of Channels"<<endl;
        continue;
      }
      
      //Make sure the event is not saturated
      if (event->event_info.saturated == true)
        continue;
        
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
        
      Double_t total_s1_corr = total_s1*ds50analysis::s1_corr_factor(t_drift_max, t_drift);
      Double_t total_s1_90 = event->pulses[s1_pulse_id].param.f90*event->pulses[s1_pulse_id].param.npe;
      Double_t total_s1_20 = event->pulses[s1_pulse_id].param.f_param[1]*event->pulses[s1_pulse_id].param.npe;
      Double_t total_f90 = total_s1_90/total_s1;
      Double_t total_f20 = total_s1_20/total_s1;
        
      Double_t total_s2 = event->pulses[s2_pulse_id].param.fixed_int2;
      Double_t total_s2_corr = total_s2/TMath::Exp(-t_drift/electron_lifetime);
        
      Double_t s2_over_s1 = total_s2/total_s1;
      Double_t s2_over_s1_corr = total_s2_corr/total_s1_corr;

      Int_t max_s1_chan = -1, max_s2_chan = -1;
      Double_t max_s1 = -1., max_s2 = -1.;
      ds50analysis::max_s1_s2(event, s1_pulse_id, s2_pulse_id, max_s1_chan, max_s1, max_s2_chan, max_s2);

      Double_t max_s1_frac = max_s1/total_s1;
      Double_t max_s2_frac = max_s2/total_s2;
      Double_t max_s1_peak_time =  event->getChannelByID(max_s1_chan).pulses[s1_pulse_id].param.peak_time;
      
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
        
        
        
        
      if (total_f90 > (0.45 + 0.5*TMath::Exp(-total_s1/50)) && total_s1 > 40 && total_s1 < 6000)
        {
          outfile<<"Run: "<<event->event_info.run_id<<" Event: "<<event->event_info.event_id<<std::endl;
          outfile<<"S1: "<<total_s1<<" S2: "<<total_s2<<std::endl;
          outfile<<"F20: "<<total_f20<<" F90: "<<total_f90<<" S2/S1: "<<s2_over_s1<<" S2/S1 Corr: "<<s2_over_s1_corr<<std::endl;
          outfile<<"Max S1 fraction: "<<max_s1_frac <<" Max S1 Channel: "<<max_s1_chan<<std::endl;
          outfile<<"Max S2 fraction: "<<max_s2_frac <<" Max S2 Channel: "<<max_s2_chan<<std::endl;
          outfile<<"TDrift: "<<t_drift<<std::endl<<std::endl;
          outfile<<"xpos: "<<xpos<<std::endl<<std::endl;
          outfile<<"ypos: "<<ypos<<std::endl<<std::endl;
        }



      // Cut on large max_s1_frac. Threshold is defined bin by bin.
      if (ds50analysis::large_max_s1_frac(total_s1, t_drift, max_s1))
        continue;

      ///////////////////////
      //  FILL HISTOGRAMS  //
      ///////////////////////
      
      for (int j = 0; j < n_hists; j++)
        {
          if (total_s1_corr > j*5 && total_s1_corr<500.)
            logs2s1_corr_f90_hist[j]->Fill(total_f90, TMath::Log10(s2_over_s1_corr));
          if (total_s1_corr >= j*5 && total_s1_corr < (j+1)*5){
            max_frac_t_drift_hist[j]->Fill(t_drift, max_s1_frac);
            if(total_s1>60. && total_s1<70. && t_drift>70. && t_drift<250 && max_s1_frac>0.2){
              cout << "runid="<< event->event_info.run_id
                   <<" event=" <<  event->event_info.event_id
                   <<" frac="<< max_s1_frac << endl;
            }
          }
        }
      max_frac_s1_hist->Fill(total_s1, max_s1_frac);

      
      
      total_s1_hist               ->Fill(total_s1);
      total_s1_corr_hist          ->Fill(total_s1_corr);
      total_s2_hist               ->Fill(total_s2);
      total_s2_corr_hist          ->Fill(total_s2_corr);
      total_f90_hist              ->Fill(total_f90);
      total_s1_f90_hist           ->Fill(total_s1, total_f90);
      total_s1_corr_f90_hist      ->Fill(total_s1_corr, total_f90);
      total_s1_corr_f20_hist      ->Fill(total_s1_corr, total_f20);
      max_frac_max_peak_time_hist ->Fill(max_s1_frac, max_s1_peak_time);
      s1_max_peak_time_hist       ->Fill(max_s1, max_s1_peak_time);
      t_drift_s1_corr_hist        ->Fill(total_s1_corr, t_drift);
      logs2s1_corr_s1_corr_hist   ->Fill(total_s1_corr, TMath::Log10(s2_over_s1_corr));
      s2_s1_corr_hist             ->Fill(total_s1_corr, total_s2_corr);
      max_s2_chan_hist            ->Fill(max_s2_chan);
    }//End loop over events
    
  std::cout<<"Run time: "<<LivetimeTotal+InhibitTimeTotal
           <<" Live time: "<<LivetimeTotal
           <<" Inhibit time: "<<InhibitTimeTotal
           <<std::endl;
  total_livetime->SetBinContent(1,LivetimeTotal);
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
  s2_s1_corr_hist->Write();
  total_f90_hist->Write();
  total_s1_f90_hist->Write();
  total_s1_corr_f20_hist->Write();
  total_s1_corr_f90_hist->Write();
  t_drift_s1_corr_hist->Write();
  max_frac_s1_hist->Write();
  s1_max_peak_time_hist->Write();
  max_frac_max_peak_time_hist->Write();
  logs2s1_corr_s1_corr_hist->Write();
  total_livetime->Write();
  max_s2_chan_hist->Write();
  for (int j = 0; j < n_hists; j++)
    {
      logs2s1_corr_f90_hist[j]->Write();
      max_frac_t_drift_hist[j]->Write();
    }
    
    
  //new TCanvas();
  //total_s1_hist->Draw("");
  //new TCanvas();
  //total_s1_corr_hist->Draw("");
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

