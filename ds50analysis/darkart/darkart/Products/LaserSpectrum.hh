/*

  Definition of the LaserSpectrum class

  AFan 2013-07-01

*/


#ifndef darkart_Products_LaserSpectrum_hh
#define darkart_Products_LaserSpectrum_hh

#include <vector>
#include "TH1F.h"
#include "TBox.h"
#include "TLine.h"
#include "TGraph.h"
#include "TFitResultPtr.h"
#include "Rtypes.h"
#include "cetlib/map_vector.h"


// Forward declarations
namespace darkart 
{
    class LaserSpectrum;
    typedef cet::map_vector<LaserSpectrum*> LaserSpectrumMap;
    class LaserResults;
    struct LaserResultsData;
}

struct darkart::LaserResultsData
{
    LaserResultsData():
	run(-9999), channel(-9999),
	int_shift_correction(-9999), int_shift_correction_unc2(-9999),
	amp_shift_correction(-9999), amp_shift_correction_unc2(-9999),
	amplitude_threshold(-9999), laser_amp_frac(-9999),
	pedestal_amp_corr(-9999),laser_amp_frac_corr(-9999),
	occupancy(-9999), occupancy_unc2(-9999),
	pedestal_mean(-9999), pedestal_mean_unc2(-9999), pedestal_var(-9999),
	total_mean(-9999), total_mean_unc2(-9999), total_var(-9999),
	gauss_center(-9999), gauss_sigma(-9999),
	ser_mean(-9999), ser_mean_unc2(-9999), ser_sigma(-9999)
	{}
    
    int run;
    int channel;
    double int_shift_correction;
    double int_shift_correction_unc2;
    double amp_shift_correction;
    double amp_shift_correction_unc2;

    double amplitude_threshold;
    double laser_amp_frac;
    double pedestal_amp_corr;
    double laser_amp_frac_corr;
    
    double occupancy;
    double occupancy_unc2;
    double pedestal_mean;
    double pedestal_mean_unc2;
    double pedestal_var;
    double total_mean;
    double total_mean_unc2;
    double total_var;
  
    double gauss_center;
    double gauss_sigma;

    double ser_mean;
    double ser_mean_unc2;
    double ser_sigma;

};

class darkart::LaserResults
{
public:
    
    LaserResults(): data(),
		    ser_mean_line(), amp_thresh_line(),
		    laser_int_sub_hist(), laser_amp_sub_hist(),
		    ped_int_scaled_hist(), ped_amp_scaled_hist()
	{
	    int_shift_chi2_graph.SetName("int_shift_chi2_graph");
	    int_shift_chi2_graph.SetMarkerStyle(20);
	    amp_shift_chi2_graph.SetName("amp_shift_chi2_graph");
	    amp_shift_chi2_graph.SetMarkerStyle(20);
	}
    
    LaserResultsData data;
    TLine ser_mean_line;
    TLine amp_thresh_line;

    TH1F laser_int_sub_hist, laser_amp_sub_hist;
    TH1F ped_int_scaled_hist, ped_amp_scaled_hist;
    TGraph int_shift_chi2_graph, amp_shift_chi2_graph;

    void DrawSubtractedIntegral()
    { 
	laser_int_sub_hist.Rebin(32);
	laser_int_sub_hist.Draw("");
	ser_mean_line.Draw("SAME");
    }
    void DrawSubtractedAmplitude()
    { 
	laser_amp_sub_hist.Rebin(80);
	laser_amp_sub_hist.Draw("");
	amp_thresh_line.Draw("SAME");
    }
  
};

class darkart::LaserSpectrum
{
    
public:
    
    TH1F laser_int_hist, laser_amp_hist,
	pedestal_int_hist, pedestal_amp_hist,
	threshold_hist;
    
    LaserResults results;
    
    LaserSpectrum(Int_t nbins,Double_t xlow,Double_t xup) :
	laser_int_hist("laser_int_hist","",nbins,xlow,xup), 
	laser_amp_hist("laser_amp_hist","",104000,-10,250),
	pedestal_int_hist("pedestal_int_hist","",nbins,xlow,xup), //NOTE: Pedestal histograms assumed to have exact same binning as laser couterparts
	pedestal_amp_hist("pedestal_amp_hist","",104000,-10,250),
	threshold_hist("threshold_hist","",nbins,xlow,xup),
	results()
    {
	laser_int_hist.Sumw2(true);
	laser_amp_hist.Sumw2(true);
	pedestal_int_hist.Sumw2(true);
	pedestal_amp_hist.Sumw2(true);
	
	laser_int_hist.SetStats(true);
	laser_amp_hist.SetStats(true);
	pedestal_int_hist.SetStats(true);
	pedestal_amp_hist.SetStats(true);
	
	laser_int_hist.StatOverflows();
	laser_amp_hist.StatOverflows();
	pedestal_int_hist.StatOverflows();
	pedestal_amp_hist.StatOverflows();
	
	laser_int_hist.SetLineColor(kBlack);
	laser_amp_hist.SetLineColor(kBlack);
	pedestal_int_hist.SetLineColor(kRed);
	pedestal_amp_hist.SetLineColor(kRed);
	threshold_hist.SetLineStyle(7);
    }
    LaserSpectrum():laser_int_hist(), laser_amp_hist(),
		    pedestal_int_hist(), pedestal_amp_hist(),
		    threshold_hist(), results()
    {}
  
    virtual ~LaserSpectrum(){}
    void Draw(bool draw_result)
    {   
	//laser_int_hist.Rebin(32);
	laser_int_hist.Draw("HIST");
	laser_int_hist.Draw("FUNC SAME");
	if (draw_result)
	    results.ser_mean_line.Draw("SAME");
    }
    void DrawIntegrals()
    { 
	//laser_int_hist.Rebin(32);
	laser_int_hist.Draw("HIST");
	pedestal_int_hist.Rebin(32);
	pedestal_int_hist.Draw("HIST SAMES");
    }
    void DrawAmplitudes()
    {   
	laser_amp_hist.Rebin(32);
	laser_amp_hist.Draw("HIST");
	pedestal_amp_hist.Rebin(32);
	pedestal_amp_hist.Draw("HIST SAMES");
    }
    
    
};




#endif
