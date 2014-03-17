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
#include "merge_lib.hh"

using namespace std;
using namespace darkart;

// Forward declaration
void LoopOverChain(TChain* tpc_chain, TChain* od_chain, TString outFileName = "analysis.root");
bool LSVMultiplicityCut(float charge, float height, float multiplicity){
  //  return height/multiplicity < 2.525e7 + sqrt(1.475e14 + 2.075e12*(charge - 14.473)*(charge - 14.473));
  return height/multiplicity < 2.563e7 + sqrt(1.574e14 + 1.39e12*(charge - 14.4)*(charge - 14.4));
}

Double_t GetPromptLSVCharge(int nclusters, float *charge, float *start, float *height, float *multiplicity){
  // Returns the charge in a prompt window around a TPC event
  float offset = 3778.; // ns
  float promptWindow = 32.; // ns
  float promptCharge = 0;
  for(int clust_num = 0; clust_num < nclusters; clust_num++){
    if(start[clust_num]-offset > promptWindow/2.)
      break;
    if(start[clust_num]-offset < -1*promptWindow/2.)
      continue;

    if(LSVMultiplicityCut(charge[clust_num], height[clust_num], multiplicity[clust_num]))
      promptCharge += charge[clust_num];
  }

  return (Double_t)promptCharge;
}

Double_t GetMaxLSVClusterCharge(int nclusters, float *charge){
  // Returns the charge of the cluster with the greatest charge within the veto gate
  float maxCharge = 0;
  for(int clust_num = 0; clust_num < nclusters; clust_num++){
    if(charge[clust_num] > maxCharge)
      maxCharge = charge[clust_num];
  }

  return (Double_t)maxCharge;
}

void analysis(TString Inputfilelist, TString outFileName) {

  TChain* tpc_chain = new TChain("treeBuilder/Events");
  TChain* od_chain = new TChain("odtree");

  std::cout<<"WARNING: Database access disabled ! Make sure to check run list manually on the ELOG"<<std::endl;

  Bool_t IsChained = AddFile2Chain(Inputfilelist, *tpc_chain, *od_chain);
  if (IsChained) LoopOverChain(tpc_chain, od_chain, outFileName);
  else std::cout<<"Cannot chain files!! Please check input file list."<<std::endl;
}

void LoopOverChain(TChain* tpc_chain, TChain *od_chain, TString outFileName)
{
  const Double_t t_drift_min = 10.; //[us]
  const Double_t t_drift_max = 373.3; //[us]
  const Double_t t_drift_delta = 10; //[us]
  const Double_t t_s3_sep_min = 372.; //372.; //[us]
  const Double_t t_s3_sep_max = 400.; //390.; //[us]
  const Double_t electron_lifetime = 4733; //3338; //3076; (old value) //[us]
  const Int_t N_CHANNELS = 38;

  // TPC Variables    
  Int_t tpc_events = tpc_chain->GetEntries();
  if (tpc_events == 0)
    {
      std::cout<<"No events loaded in TPC TChain. Aborting..."<<std::endl;
      return;
    }
  std::cout << "Total number of TPC events: "<<tpc_events << endl;
    
  EventData* event = NULL;
  tpc_chain->SetBranchAddress("EventData", &event);
    
  // OD Variables
  Int_t od_events = od_chain->GetEntries();
  if(od_events == 0)
    {
      std::cout << "No events loaded in OD TChain. Aborting..." << std::endl;
      return;
    }
  std::cout << "Total number of OD events: " << od_events << endl;

  double max_dt = 64*20.e-9; // s
  Int_t od_run=0;
  Int_t od_event=0;
  Int_t lsv_nclusters=0;
  const int od_max_nclusters = 40;
  Float_t lsv_cluster_charge[od_max_nclusters];
  Float_t lsv_cluster_start[od_max_nclusters];
  Float_t wt_total_charge = 0;
  Float_t lsv_cluster_height[od_max_nclusters];
  Float_t lsv_cluster_multiplicity[od_max_nclusters];
  UInt_t od_gps_fine_time_counter = 0;
  UShort_t od_pps_counter = 0;
  const Double_t max_lsv_charge_thresh = 150;
  const Double_t prompt_lsv_charge_thresh = 3;
  const Double_t wt_total_charge_thresh = 200;

  od_chain->SetBranchAddress("run", &od_run);
  od_chain->SetBranchAddress("event_number", &od_event);
  od_chain->SetBranchAddress("lsv_n_clusters", &lsv_nclusters);
  od_chain->SetBranchAddress("lsv_cluster_fixed_width_charge", lsv_cluster_charge);
  od_chain->SetBranchAddress("lsv_cluster_start_ns", lsv_cluster_start);
  od_chain->SetBranchAddress("wt_total_spe_charge", &wt_total_charge);
  od_chain->SetBranchAddress("lsv_cluster_height", lsv_cluster_height);
  od_chain->SetBranchAddress("lsv_cluster_max_multiplicity", lsv_cluster_multiplicity);
  od_chain->SetBranchAddress("gps_fine_time_counter", &od_gps_fine_time_counter);
  od_chain->SetBranchAddress("pps_counter", &od_pps_counter);

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
  TH1F* t_drift_hist                = new TH1F("t_drift_hist", "Drift Time; t_drift [#mus]",
                                               400, 0., 400.);
  TH1F* s1_startime_hist            = new TH1F("s1_startime_hist", "S1 Start Time; start time [#mus]",
                                               500, -1., 1.);
  TH1F* total_s1_hist               = new TH1F("total_s1_hist", "S1 Spectrum; S1 [PE]",
                                               10000, 0, 10000);
  TH1F* total_s1_corr_hist          = new TH1F("total_s1_corr_hist", "S1 Spectrum (corrected for z-dependence); total_s1_corr [PE]",
                                               10000, 0, 10000);
  TH1F* total_s2_hist               = new TH1F("total_s2_hist", "S2 Spectrum; S2 [PE]",
                                               5000, 0, 500000);
  TH1F* total_s2_corr_hist          = new TH1F("total_s2_corr_hist", "S2 Spectrum (corrected for z-dependence); total_s2_corr [PE]",
                                               5000, 0, 500000);
  TH1F* total_s1_90_hist            = new TH1F("total_s1_90_hist", "S1_{90} Spectrum; S1_{90} [PE]",
                                               2000, 0, 2000);
  TH1F* total_s1_90_corr_hist       = new TH1F("total_s1_90_corr_hist", "S1_{90} Spectrum (corrected for z-dependence); S1_{90} (corr) [PE]",
                                               2000, 0, 2000);
  TH1F* total_f90_hist              = new TH1F("total_f90_hist", "F90 Distribution; F90",
                                               110, 0, 1.3);
  TH1F* total_f90_prompt_cut_hist   = new TH1F("total_f90_prompt_cut_hist", "F90 Distribution, with Prompt Cut; F90",
                                               110, 0, 1.3);
  TH1F* veto_prompt_hist            = new TH1F("veto_prompt_hist", "Veto Prompt Spectrum; lsv_prompt_charge [PE]",
                                               4000, 0, 4000);
  TH1F* veto_prompt_3500_10000_hist = new TH1F("veto_prompt_3500_10000_hist", "Veto Prompt Spectrum for TPC Events between 3500 and 10000 PE; lsv_prompt_charge [PE]",
					       4000, 0, 4000);
  TH1F* veto_prompt_1200_2500_hist = new TH1F("veto_prompt_1200_2500_hist", "Veto Prompt Spectrum for TPC Events between 1200 and 2500 PE; lsv_prompt_charge [PE]",
					       2000, 0, 2000);
  TH1F* veto_prompt_0_1200_hist     = new TH1F("veto_prompt_0_1200_hist", "Veto Prompt Spectrum for TPC Events between 0 and 1200 PE; lsv_prompt_charge [PE]",
					       2000, 0, 2000);
  TH2F* veto_prompt_tpc_s1_hist     = new TH2F("veto_prompt_tpc_s1_hist", "Veto Prompt Charge vs. TPC S1;total_s1_corr [PE]; lsv_prompt_charge [PE]",
					       10000, 0, 10000, 2000, 0, 2000);
  TH1F* veto_max_charge_hist        = new TH1F("veto_max_charge_hist", "Veto Maximum Cluster Charge Spectrum; Charge [p.e.]",
                                               4000,0,4000);
  TH1F* total_f90_veto_cut_hist     = new TH1F("total_f90_veto_cut_hist", "F90 Distribution, with Veto Cut; F90",
                                               110, 0, 1.3);
  TH2F* total_s1_f90_hist           = new TH2F("total_s1_f90_hist", "F90 vs S1; S1 [PE]; F90",
                                               10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_hist      = new TH2F("total_s1_corr_f90_hist", "F90 vs S1 (corrected for z-dependence); total_s1_corr [PE]; F90",
                                               10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f20_hist      = new TH2F("total_s1_corr_f20_hist", "F20 vs S1 (corrected for z-dependence); total_s1_corr [PE]; F20",
                                               10000, 0, 10000, 130, 0, 1.3);
  TH2F* max_frac_s1_hist            = new TH2F("max_frac_s1_hist", "Maximum S1 fraction vs S1; S1 [PE]; max_s1/total_s1",
                                               10000, 0, 10000, 100, 0, 1.0);
  TH2F* max_frac_max_peak_time_hist = new TH2F("max_frac_max_peak_time_hist", "Max Peak Arrival Time vs Maximum S1 fraction vs ",
                                               10000, 0, 1.0, 800, 0, 8.0);
  TH2F* s1_max_peak_time_hist       = new TH2F("s1_max_peak_time_hist", "Max Peak Arrival Time vs Maximum S1",
                                               10000, 0, 10000, 100, 0, 8.0);
  TH2F* t_drift_total_s1_hist        = new TH2F("t_drift_total_s1_hist", "Drift time vs S1; total_s1 [PE]; t_drift [#mus]",
                                               6000, 0, 6000, 400, 0., 400.);
  TH2F* t_drift_s1_corr_hist        = new TH2F("t_drift_s1_corr_hist", "Drift time vs S1 (corrected for z-dependence); total_s1 [PE]; t_drift [#mus]",
                                               6000, 0, 6000, 400, 0., 400.);
  TH2F* s2_s1_corr_hist             = new TH2F("s2_s1_corr_hist", "S2 vs S1 (corrected for z-dependence); total_s1_corr [PE], total_s2_corr [PE]",
                                               6000, 0, 6000, 5000, 0, 50000);
  TH2F* logs2s1_corr_s1_corr_hist   = new TH2F("s2s1_corr_s1_corr_hist",
                                               "Log(S2/S1) vs S1, S1 (corrected for z-dependence); total_s1_corr [PE]; Log(S2/S1)",
                                               6000, 0, 6000, 200, -1, 3);

  TH2F* total_s1_f90_after_lsv_cuts_hist       = new TH2F("total_s1_f90_after_lsv_cuts_hist",
                                                          "F90 vs S1, TPC and LSV cuts; S1 [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_f90_prompt_lsv_cut_hist       = new TH2F("total_s1_f90_prompt_lsv_cut_hist",
                                                          "F90 vs S1, TPC cuts, prompt LSV cut; S1 [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_f90_max_lsv_cut_hist          = new TH2F("total_s1_f90_max_lsv_cut_hist",
                                                          "F90 vs S1, TPC cuts, max LSV cut; S1 [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_f90_after_lsv_wt_cuts_hist    = new TH2F("total_s1_f90_after_lsv_wt_cuts_hist",
                                                          "F90 vs S1, TPC cuts, TPC and LSV and WT cuts; S1 [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_after_lsv_cuts_hist  = new TH2F("total_s1_corr_f90_after_lsv_cuts_hist",
                                                          "F90 vs S1 (corrected for z-dependence), TPC and LSV cuts; total_s1_corr [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_prompt_lsv_cut_hist  = new TH2F("total_s1_corr_f90_prompt_lsv_cut_hist",
                                                          "F90 vs S1 (corrected for z-dependence), TPC cuts, prompt LSV cut; total_s1_corr [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_max_lsv_cut_hist     = new TH2F("total_s1_corr_f90_max_lsv_cut_hist",
                                                          "F90 vs S1 (corrected for z-dependence), TPC cuts, max LSV cut; total_s1_corr [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_prompt_lsv_cut_c_hist= new TH2F("total_s1_corr_f90_prompt_lsv_cut_c_hist",
                                                          "F90 vs S1 (corrected for z-dependence), TPC cuts, NOT passing prompt LSV cut; total_s1_corr [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_max_lsv_cut_c_hist   = new TH2F("total_s1_corr_f90_max_lsv_cut_c_hist",
                                                          "F90 vs S1 (corrected for z-dependence), TPC cuts, NOT passing max LSV cut; total_s1_corr [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH2F* total_s1_corr_f90_after_lsv_wt_cuts_hist=new TH2F("total_s1_corr_f90_after_lsv_wt_cuts_hist",
                                                          "F90 vs S1 (corrected for z-dependence), TPC and LSV and WT cuts; total_s1_corr [PE]; F90",
                                                          10000, 0, 10000, 130, 0, 1.3);
  TH1F* total_s1_corr_prompt_lsv_cut_hist      = new TH1F("total_s1_corr_prompt_lsv_cut_hist",
                                                          "S1 (corrected for z-dependence), TPC cuts, prompt LSV cut; total_s1_corr [PE]",
                                                          10000, 0, 10000);
  TH1F* total_s1_corr_prompt_lsv_cut_c_hist    = new TH1F("total_s1_corr_prompt_lsv_cut_c_hist",
                                                          "S1 (corrected for z-dependence), TPC cuts, NOT passing prompt LSV cut; total_s1_corr [PE]",
                                                          10000, 0, 10000);

  
  TH1F* total_livetime              = new TH1F("total_livetime", "total_livetime", 1, 0, 1);
  TH1F* total_livetime_after_veto_cuts = new TH1F("total_livetime_after_veto_cuts", "total_livetime_after_veto_cuts", 1, 0, 1);
  TH1F* total_tpc_events            = new TH1F("total_tpc_events", "Total TPC events in chain", 1,0,1);
  TH1F* total_od_events             = new TH1F("total_od_events", "Total OD events in chain", 1,0,1);
  TH1F* total_tpc_events_used       = new TH1F("total_tpc_events_used", "Total TPC events used", 1,0,1);
  TH1F* total_od_events_used        = new TH1F("total_od_events_used", "Total OD events used", 1,0,1);
  
  const int n_hists = 70;
  const int s1_slice_width = 10;
  TH2F** logs2s1_corr_f90_hist = new TH2F*[n_hists];
  TH2F** max_frac_t_drift_hist = new TH2F*[n_hists];
  for (int i = 0; i < n_hists; i++)
    {
      std::ostringstream name, title;
      name<<"logs2s1_f90_hist_"<<i;
      title<<"Log(S2/S1) (both corr.) vs F90 for "<<i*s1_slice_width<<" < total_s1_corr [p.e.] < "<<(i+1)*s1_slice_width<<"; F90; Log(S2/S1)";
      logs2s1_corr_f90_hist[i] = new TH2F(name.str().c_str(), title.str().c_str(), 120, 0.0, 1.2, 200, -1, 3);

      name.str("");title.str("");
      name<<"max_frac_t_drift_hist_"<<i;
      title<<"Max S1 Fraction vs Drift time for "<<(i)*s1_slice_width<<" < S1 [p.e.] < "<<(i+1)*s1_slice_width;
      max_frac_t_drift_hist[i] = new TH2F(name.str().c_str(), title.str().c_str(), 400, 0, 400, 100, 0, 1.0);
    }
    
  Double_t LivetimeTotal(0.), InhibitTimeTotal(0.), LivetimeAfterVetoCuts(0.), NumODEventsUsed(0.), NumTPCEventsUsed(0.);

  //////////////////////////////////////////////////////////////////////////
  /////////////////     BEGIN LOOP OVER EVENTS       ///////////////////////
  //////////////////////////////////////////////////////////////////////////
  //tpc_events = tpc_events - 50; // Skip last few events because end of some (very few) runs are problematic.
  int m = 0; // m is the event index into the OD tree. Note that TPC event numbers are 1-indexed, while OD event numbers are 0-indexed
  int old_tpc_run = 0, old_od_run = 0;
  for (int n = 0; n < tpc_events; n++)
    {
      if(m > od_events-1){
	std::cout << "Ran out of OD events at OD event # " << m << "  and TPC event number " << n << std::endl;
	break; // We've run out of OD events
      }
        
      //Load the event
      tpc_chain->GetEntry(n);
      od_chain->GetEntry(m);
    
      // Check to see if you've just opened a new TPC data file
      if(event->event_info.run_id != old_tpc_run)
	{
	  std::cout << "New TPC Run found, start new OD Run, too: starting new set of runs" << std::endl;
	  //If you didn't also just open a new OD data file, burn through the chain until you reach the next one
	  while(old_od_run == od_run)
	    od_chain->GetEntry(++m);

	  old_od_run = od_run;
	  old_tpc_run = event->event_info.run_id;
	  std::cout << "\tOpening TPC Run " << old_tpc_run << " and OD Run " << old_od_run << std::endl;
	}
      // Check to see if you've just opened a new OD data file
      else if(od_run != old_od_run)
	{
	  std::cout << "New OD Run found, start new TPC Run, too: starting new set of runs" << std::endl;
	  //If you didn't also just open a new OD data file, burn through the chain until you reach the next one
	  while(old_tpc_run == event->event_info.run_id)
	    tpc_chain->GetEntry(++n);

	  old_od_run = od_run;
	  old_tpc_run = event->event_info.run_id;
	  std::cout << "\tOpening TPC Run " << old_tpc_run << " and OD Run " << old_od_run << std::endl;
	}
        
      Int_t run_id = event->event_info.run_id;
      Double_t acquiWindow = (1.*event->sumchannel.channel.nsamps)/(1.*event->sumchannel.channel.sample_rate);
      if (n % 10000 == 0)
	{
	  std::cout<<"Processing\t TPC Event: "<<n<<"/"<<tpc_events<<std::endl;


      if(n==0)std::cout<<"acquisition window [us]: "<<acquiWindow<<std::endl;

//      if(acquiWindow < 565 || acquiWindow > 575)std::cout<<"Livetime cut does not work with this acquisition window!!!!"<<std::endl;
	}

      
      // There is some mismatch between events, so check that the events in both trees agree
      int test_m = m;
      int start_od_run = od_run;
      double gps_diff = ((double)od_pps_counter - (double)event->event_info.gps_coarse +
                         ((double)od_gps_fine_time_counter - (double)event->event_info.gps_fine)*20.e-9);
      
      if (gps_diff >= max_dt || gps_diff <= 0) {
        while (gps_diff >= max_dt || gps_diff <= 0) {
	  test_m++;
	  if(test_m >= od_events || od_run != start_od_run){
	    std::cout << "Reach end of run: " << "start_od_run = " << start_od_run << " and od_run = " << od_run << std::endl;
	    if(od_run != start_od_run)
	      m = test_m-1;
	    break;
	  }
	  od_chain->GetEntry(test_m);
          gps_diff = ((double)od_pps_counter - (double)event->event_info.gps_coarse +
                      ((double)od_gps_fine_time_counter - (double)event->event_info.gps_fine)*20.e-9);
	}
        if (test_m >= od_events && gps_diff < 0) {
          std::cout << "Current TPC event: run "<<event->event_info.run_id<<", event "<<event->event_info.event_id<<", gps time "
                    << event->event_info.gps_coarse+event->event_info.gps_fine*20.e-9<<"\n"
                    << "Final OD event in chain: run "<<od_run<<", event "<<od_event<<", gps time "
                    << od_pps_counter+od_gps_fine_time_counter*20.e-9
                    << ". No more events will match. Exiting loop."<< std::endl;
          break;
        }
	if(test_m >= od_events || od_run != start_od_run){
	  {
            std::cout << "Could not find an OD event to match TPC event " << event->event_info.event_id
                      << " in run " << event->event_info.run_id << " ... skipping event" << std::endl;

            
	    continue;
	  }
	}

      }
      m = ++test_m;
      
      if (n % 10000 == 0)
	{
	  std::cout<<"\t\t OD Event : " << m << "/" << od_events << std::endl;
	}

      NumODEventsUsed++;
      NumTPCEventsUsed++;
      if(event->event_info.live_time_20ns*20.*1.e-9 < 1.)
        { //long Livetime event cut
          LivetimeTotal+=event->event_info.live_time_20ns*20.*1.e-9; // in second
          InhibitTimeTotal+=event->event_info.incremental_inhibit_time_20ns*20.*1.e-9; // in second
        }// This should be before any event cuts!!
      else continue; // not include long livetime runs



      /////////////////////////
      //  APPLY BASIC CUTS   //
      /////////////////////////

      //CALCULATE OD PARAMETERS
      
      // The total LSV cluster charge arriving in a 20ns prompt window
      Double_t prompt_lsv_charge = GetPromptLSVCharge(lsv_nclusters, lsv_cluster_charge, lsv_cluster_start, lsv_cluster_height, lsv_cluster_multiplicity);
      // The maximum charge found in any LSV clusters in the 45 us veto gate
      Double_t max_lsv_charge = GetMaxLSVClusterCharge(lsv_nclusters, lsv_cluster_charge); 

      
      if ( !(max_lsv_charge > max_lsv_charge_thresh || prompt_lsv_charge > prompt_lsv_charge_thresh || wt_total_charge > wt_total_charge_thresh ) ) {
        if(event->event_info.live_time_20ns*20.*1.e-9 < 1.)
          { //long Livetime event cut
            LivetimeAfterVetoCuts+=event->event_info.live_time_20ns*20.*1.e-9; // in second
          }// This should be before any event cuts!!
        else continue; // not include long livetime runs
      }

      //if ( wt_total_charge > 0 ) continue;

      // Check for expected number of channels
      if ((int)event->channels.size() != N_CHANNELS){
        cout << "Event=" << event->event_info.event_id<<" has LOWER # of Channels"<<endl;
        continue;
      }

      //Make sure the baseline was found on the sum channel
      if (event->sumchannel.baseline.found_baseline == false)
        continue;
        
      //PULSE IDENTIFICATION
      Int_t n_phys_pulses = -1, s1_pulse_id = -1, s2_pulse_id = -1;
      ds50analysis::identify_pulses(event, n_phys_pulses, s1_pulse_id, s2_pulse_id, t_s3_sep_min, t_s3_sep_max);
      
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
      Double_t total_s1_90_corr = total_s1_90*ds50analysis::s1_corr_factor(t_drift_max, t_drift);
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

      //Remove events triggered on S3
      Double_t dt_usec = (event->event_info.live_time_20ns + event->event_info.incremental_inhibit_time_20ns)*20.*1.e-3; //time from previous trigger
      Double_t inhibit_window_us = (acquiWindow > 150.)? (acquiWindow - 60.)*2. + 50. : 0.; // in us // 60 for large S2 signal after max drift time. 50 for buffer for S3
      if (dt_usec < inhibit_window_us)
        continue;

#define S1SATURATION
#ifndef S1SATURATION
     //Make sure the event is not saturated
      if (event->event_info.saturated == true)
        continue;
#else
      //Make sure the s1 is not saturated
      if (event->pulses[s1_pulse_id].param.peak_saturated == true)
        continue;
#endif

      //Make sure the S1 pulse is where you expect in the trigger window
      if (s1_start_time < -0.25 || s1_start_time > -0.18 )
        continue;

      t_drift_hist->Fill(t_drift);
      
      //Remove very slow events (triggered on S2)
      if (total_f90 < 0.1)
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
     

      // Apply veto cuts and fill subset of histos separately
      if (max_lsv_charge <= max_lsv_charge_thresh) {
        total_s1_f90_max_lsv_cut_hist->Fill(total_s1, total_f90);
        total_s1_corr_f90_max_lsv_cut_hist->Fill(total_s1_corr, total_f90);
      }
      else {
        total_s1_corr_f90_max_lsv_cut_c_hist->Fill(total_s1_corr, total_f90);
      }
      if (prompt_lsv_charge <= prompt_lsv_charge_thresh) {
        total_s1_f90_prompt_lsv_cut_hist->Fill(total_s1, total_f90);
        total_s1_corr_f90_prompt_lsv_cut_hist->Fill(total_s1_corr, total_f90);
        total_s1_corr_prompt_lsv_cut_hist->Fill(total_s1_corr);
      }
      else {
        total_s1_corr_f90_prompt_lsv_cut_c_hist->Fill(total_s1_corr, total_f90);
        total_s1_corr_prompt_lsv_cut_c_hist->Fill(total_s1_corr);
      }


      if (max_lsv_charge <= max_lsv_charge_thresh && prompt_lsv_charge <= prompt_lsv_charge_thresh) {
        total_s1_f90_after_lsv_cuts_hist->Fill(total_s1, total_f90);
        total_s1_corr_f90_after_lsv_cuts_hist->Fill(total_s1_corr, total_f90);
        if ( wt_total_charge <= wt_total_charge_thresh ) {
          total_s1_f90_after_lsv_wt_cuts_hist->Fill(total_s1, total_f90);
          total_s1_corr_f90_after_lsv_wt_cuts_hist->Fill(total_s1_corr, total_f90);
        }
      }
      
      ///////////////////////
      //  FILL HISTOGRAMS  //
      ///////////////////////
      
      for (int j = 0; j < n_hists; j++)
        {
          //if (total_s1_corr > j*5 && total_s1_corr<500.)
          if (total_s1_corr > j*s1_slice_width && total_s1_corr <= (j+1)*s1_slice_width) {
            logs2s1_corr_f90_hist[j]->Fill(total_f90, TMath::Log10(s2_over_s1_corr));
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
      total_s1_90_hist            ->Fill(total_s1_90);
      total_s1_90_corr_hist       ->Fill(total_s1_90_corr);
      total_f90_hist              ->Fill(total_f90);
      total_s1_f90_hist           ->Fill(total_s1, total_f90);
      total_s1_corr_f90_hist      ->Fill(total_s1_corr, total_f90);
      total_s1_corr_f20_hist      ->Fill(total_s1_corr, total_f20);
      max_frac_max_peak_time_hist ->Fill(max_s1_frac, max_s1_peak_time);
      s1_max_peak_time_hist       ->Fill(max_s1, max_s1_peak_time);
      t_drift_total_s1_hist       ->Fill(total_s1, t_drift);
      t_drift_s1_corr_hist        ->Fill(total_s1_corr, t_drift);
      logs2s1_corr_s1_corr_hist   ->Fill(total_s1_corr, TMath::Log10(s2_over_s1_corr));
      s2_s1_corr_hist             ->Fill(total_s1_corr, total_s2_corr);
      veto_prompt_hist            ->Fill(prompt_lsv_charge);
      if(prompt_lsv_charge > 0){
	if(total_s1_corr > 0 && total_s1_corr < 1200)
	  veto_prompt_0_1200_hist->Fill(prompt_lsv_charge);
	if(total_s1_corr > 1200 && total_s1_corr < 2500)
	  veto_prompt_1200_2500_hist->Fill(prompt_lsv_charge);
	if(total_s1_corr > 3500 && total_s1_corr < 10000)
	  veto_prompt_3500_10000_hist->Fill(prompt_lsv_charge);
	veto_prompt_tpc_s1_hist->Fill(total_s1_corr, prompt_lsv_charge);
      }
      veto_max_charge_hist        ->Fill(max_lsv_charge);
      if(prompt_lsv_charge > prompt_lsv_charge_thresh) 
	total_f90_prompt_cut_hist ->Fill(total_f90);
      if(max_lsv_charge > max_lsv_charge_thresh)
	total_f90_veto_cut_hist   ->Fill(total_f90);

    }//End loop over events
    
  std::cout<<"Run time: "<<LivetimeTotal+InhibitTimeTotal << endl
           <<" Live time: "<<LivetimeTotal << endl
	   <<" Live time after vetoing on OD " << LivetimeAfterVetoCuts << endl
           <<" Inhibit time: "<<InhibitTimeTotal << endl
	   <<" Lost " << (1.- NumODEventsUsed/od_events)*100. << " of all veto events" << endl
	   <<" Lost " << (1.- NumTPCEventsUsed/tpc_events)*100. << " of all tpc events"<< endl
           <<std::endl;
  total_livetime->SetBinContent(1, LivetimeTotal);
  total_livetime_after_veto_cuts->SetBinContent(1, LivetimeAfterVetoCuts);
  total_tpc_events->SetBinContent(1, tpc_events);
  total_od_events->SetBinContent(1, od_events);
  total_tpc_events_used->SetBinContent(1, NumTPCEventsUsed);
  total_od_events_used->SetBinContent(1, NumODEventsUsed);
  
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
  total_s1_90_hist->Write();
  total_s1_90_corr_hist->Write();
  s2_s1_corr_hist->Write();
  total_f90_hist->Write();
  total_f90_prompt_cut_hist->Write();
  total_f90_veto_cut_hist->Write();
  veto_prompt_hist->Write();
  veto_prompt_0_1200_hist->Write();
  veto_prompt_1200_2500_hist->Write();
  veto_prompt_3500_10000_hist->Write();
  veto_prompt_tpc_s1_hist->Write();
  veto_max_charge_hist->Write();
  total_s1_f90_hist->Write();
  total_s1_corr_f20_hist->Write();
  total_s1_corr_f90_hist->Write();
  t_drift_total_s1_hist->Write();
  t_drift_s1_corr_hist->Write();
  max_frac_s1_hist->Write();
  s1_max_peak_time_hist->Write();
  max_frac_max_peak_time_hist->Write();
  total_livetime->Write();
  total_livetime_after_veto_cuts->Write();
  logs2s1_corr_s1_corr_hist->Write();
  total_s1_f90_after_lsv_cuts_hist->Write();
  total_s1_f90_prompt_lsv_cut_hist->Write();
  total_s1_f90_max_lsv_cut_hist->Write();
  total_s1_f90_after_lsv_wt_cuts_hist->Write();
  total_s1_corr_f90_after_lsv_cuts_hist->Write();
  total_s1_corr_f90_prompt_lsv_cut_hist->Write();
  total_s1_corr_f90_max_lsv_cut_hist->Write();
  total_s1_corr_f90_prompt_lsv_cut_c_hist->Write();
  total_s1_corr_f90_max_lsv_cut_c_hist->Write();
  total_s1_corr_f90_after_lsv_wt_cuts_hist->Write();
  total_s1_corr_prompt_lsv_cut_hist->Write();
  total_s1_corr_prompt_lsv_cut_c_hist->Write();
  
  total_tpc_events->Write();
  total_od_events->Write();
  total_tpc_events_used->Write();
  total_od_events_used->Write();
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
  string defaultFName = "analysis.root";
  if ( argc == 2 ) {
    std::cout << "\n==========> analysis with file list <=============" << std::endl;
    std::cout << "Using default output file name : " << defaultFName << std::endl;
    analysis(argv[1], defaultFName.c_str());
  } else if ( argc == 3 ) {
    std::cout << "\n==========> analysis with file list & outputname <=============" << std::endl;
    analysis(argv[1], argv[2]);
  } else {
    std::cout << "Usage:" <<argc<< std::endl;
    std::cout << "./merge filelist.txt" << std::endl;
    std::cout << "OR ./merge filelist.txt output.root" << std::endl;
    return 0;
  }


  std::cout << "==> Application finished." << std::endl;

  return 0;
}
#endif /* __CINT __ */

