/*
  Macro to develop t_drift and total_s1 dependence of max_s1/total_s1 cut.
  Forked from stripped down version of analysis_jan2014note.C

  Alden Fan -- 2013-02-03
  
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
  //const Double_t electron_lifetime = 3330; //[us]
  const Int_t N_CHANNELS = 38;

  //const Double_t s1_corr_limit = 2000.;
  //const Int_t n_s1_slices = 20;
  //const Int_t s1_slice_width = s1_corr_limit / n_s1_slices;
  const Int_t n_tdrift_slices = 37;
  const Double_t tdrift_slice_width = 370 / n_tdrift_slices;
  

    
  Int_t tpc_events = tpc_chain->GetEntries();
  if (tpc_events == 0)
    {
      std::cout<<"No events loaded in TChain. Aborting..."<<std::endl;
      return;
    }
  std::cout << "Total number of events: "<<tpc_events << endl;
    
  EventData* event = NULL;
  tpc_chain->SetBranchAddress("EventData", &event);
  
    
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
  TH1F* total_f90_hist              = new TH1F("total_f90_hist", "F90 Distribution; F90", 110, 0, 1.3);
  TH2F* total_s1_f90_hist           = new TH2F("total_s1_f90_hist", "F90 vs S1; S1 [p.e.]; F90", 10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_hist      = new TH2F("total_s1_corr_f90_hist", "F90 vs S1 (corrected for z-dependence); S1 [p.e.]; F90",
                                               10000, 0, 10000, 130, 0, 1.3);
  
  
  
  TH2F** max_s1_frac_total_s1_hists = new TH2F*[n_tdrift_slices];
  TH1F** max_s1_frac_hists = new TH1F*[n_tdrift_slices];
  TH1F*** max_s1_frac_hists_total_s1 = new TH1F**[n_tdrift_slices*15];
  

  for (int i=0; i<n_tdrift_slices; i++) {
    std::ostringstream name, title;
    name<<"max_s1_frac_total_s1_hist_"<<i;
    title<<"max_s1/total_s1 vs total_s1, "
         <<i*tdrift_slice_width<<" < T_drift < "<<(i+1)*tdrift_slice_width<<"; total_s1 [pe]; max_s1 / total_s1";
    max_s1_frac_total_s1_hists[i] = new TH2F(name.str().c_str(), title.str().c_str(), 500, 0, 500, 400, 0, 1);

    name.str(""); title.str("");
    name<<"max_s1_frac_hist_"<<i;
    title<<"max_s1/total_s1, 200<total_s1<500, "
         <<i*tdrift_slice_width<<" < T_drift < "<<(i+1)*tdrift_slice_width<<"; max_s1 / total_s1";
    max_s1_frac_hists[i] = new TH1F(name.str().c_str(), title.str().c_str(), 800, 0, 1);

    max_s1_frac_hists_total_s1[i] = new TH1F*[15];
    for (int j=0; j<15; j++) {
      name.str(""); title.str("");
      name<<"max_s1_frac_hist_"<<i<<"_total_s1_"<<j;
      title<<"max_s1/total_s1, "<<j*10<<" < total S1 < "<<(j+1)*10<<", "
           <<i*tdrift_slice_width<<" < T_drift < "<<(i+1)*tdrift_slice_width<<"; max_s1 / total_s1";
      max_s1_frac_hists_total_s1[i][j] = new TH1F(name.str().c_str(), title.str().c_str(), 400, 0, 1);
    }
  }

  //////////////////////////////////////////////////////////////////////////
  /////////////////     BEGIN LOOP OVER EVENTS       ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  //tpc_events = tpc_events - 50; // Skip last few events because end of some (very few) runs are problematic.
  for (int n = 0; n < tpc_events; n++)
  //for (int n=0; n<20000; n++)
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
        
      Double_t total_s1_corr = total_s1*ds50analysis::s1_corr_factor(t_drift_max, t_drift);
      Double_t total_s1_90 = event->pulses[s1_pulse_id].param.f90*event->pulses[s1_pulse_id].param.npe;
      Double_t total_f90 = total_s1_90/total_s1;

      Int_t max_s1_chan = -1, max_s2_chan = -1;
      Double_t max_s1 = -1., max_s2 = -1.;
      ds50analysis::max_s1_s2(event, s1_pulse_id, s2_pulse_id, max_s1_chan, max_s1, max_s2_chan, max_s2);
      
      Double_t max_s1_frac = max_s1/total_s1;

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

      /*
      //Remove events concentrated on a single PMT
      if (max_s1/total_s1 > 0.4)
        continue;
      
      //Remove events near grid or cathode
      if (t_drift < t_drift_min || t_drift > t_drift_max - t_drift_delta)
        continue;
      
      // Cut on large max_s1_frac. Threshold is defined bin by bin.
      if (ds50analysis::large_max_s1_frac(total_s1, t_drift, max_s1))
        continue;
      */

      
      //Make sure the event is not saturated
      if (event->event_info.saturated == true)
        continue;


      
      ///////////////////////
      //  FILL HISTOGRAMS  //
      ///////////////////////
      
      
      total_s1_hist               ->Fill(total_s1);
      total_s1_corr_hist          ->Fill(total_s1_corr);
      total_f90_hist              ->Fill(total_f90);
      total_s1_f90_hist           ->Fill(total_s1, total_f90);
      total_s1_corr_f90_hist      ->Fill(total_s1_corr, total_f90);

      for (int i=0; i<n_tdrift_slices; i++) {
        if (t_drift > i*tdrift_slice_width && t_drift <= (i+1)*tdrift_slice_width) {
          max_s1_frac_total_s1_hists[i]->Fill(total_s1, max_s1_frac);
          if (total_s1 > 200 && total_s1 < 500)
            max_s1_frac_hists[i]->Fill(max_s1_frac);
          for (int j=0; j<15; j++) {
            if (total_s1 > j*10 && total_s1 <= (j+1)*10)
              max_s1_frac_hists_total_s1[i][j]->Fill(max_s1_frac);
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
  total_f90_hist->Write();
  total_s1_f90_hist->Write();
  total_s1_corr_f90_hist->Write();


  for (int i=0; i<n_tdrift_slices; i++) {
    std::ostringstream dirname;
    dirname<<"tdrift_slice_"<<i;
    f->mkdir(dirname.str().c_str());
    f->cd(dirname.str().c_str());

    max_s1_frac_total_s1_hists[i]->Write();
    max_s1_frac_hists[i]->Write();

    dirname.str("");
    dirname<<"tdrift_slice_"<<i<<"/s1_slices";
    f->mkdir(dirname.str().c_str());
    f->cd(dirname.str().c_str());
    for (int j=0; j<15; j++) {
      max_s1_frac_hists_total_s1[i][j]->Write();
    }
  }
  
  //new TCanvas();
  //total_s1_hist->Draw("");
  //new TCanvas();
  //total_s1_corr_hist->Draw("");
  f->cd();
  f->Close();

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

