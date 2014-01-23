/*

  Definition of the Pmt class

  AFan 2013-06-30

  2013-08-01 AFan -- TODO use getters to ensure that
  class cannot be instantiated in inconsistent state

 */

#ifndef darkart_Products_Pmt_hh
#define darkart_Products_Pmt_hh

#include <vector>
#include <string>
#include <iostream>
#include <cmath>


// Forward declarations
namespace darkart
{
  struct Pmt;
  typedef std::vector<Pmt> PmtVec;
}


// Holds relevant info for each PMT.
struct darkart::Pmt
{
  // Default constructor
  Pmt():
    serial_id(""),
    channel_id(-1),
    spe_mean(1.0),
    spe_sigma(0.0),
    photocathode_x(-1.0),
    photocathode_y(-1.0),
    photocathode_z(-1.0),
    //photocathode_r(-1.0),
    //photocathode_theta(0.0),
    photocathode_area(-1.0),
    qe(-1.0)
  { }

  // All the relevant info we want to store for each PMT
  std::string     serial_id; //serial number
  int            channel_id; //channelID associated with this PMT
  double           spe_mean; //mean photoelectron response read from database
  double          spe_sigma; //sigma of photoelectron response read from database
  double     photocathode_x; //x position of the center of the photocathode
  double     photocathode_y; //y position of the center of the photocathode
  double     photocathode_z; //z position of the center of the photocathode
  double     photocathode_r() const
  {
    //r (cylindrical coord) position of the center of the photocathode
    return std::sqrt(photocathode_x*photocathode_x+photocathode_y*photocathode_y);
  }
  double photocathode_theta() const
  {
    //theta (cylindrical coord) position of the center of the photocathode
    if (photocathode_x == 0 && photocathode_y == 0)
      return 0.;
    return std::atan2(photocathode_y, photocathode_x);
  }
  double  photocathode_area; //area [cm2] of the photocathode
  double                 qe; //quantum efficiency

  // AWW 24.11.13
  inline void print(int verbosity = 0) const {
    std::cout << "PMT Info:"                                                   << std::endl;
    std::cout << "----------------------------------------"                    << std::endl;
    std::cout << "Serial #:                   " << serial_id                   << std::endl;
    std::cout << "Channel ID:                 " << channel_id                  << std::endl;
    std::cout << "Quantum Efficiency:         " << qe                          << std::endl;
    std::cout << "SPE mean:                   " << spe_mean                    << std::endl;

    if(verbosity > 0){
      std::cout << "SPE sigma:                  " << spe_sigma                   << std::endl;
      std::cout << "Center X:                   " << photocathode_x              << std::endl;
      std::cout << "Center Y:                   " << photocathode_y              << std::endl;
      std::cout << "Center Z:                   " << photocathode_z              << std::endl;
      std::cout << "Center R:                   " << photocathode_r()            << std::endl;
      std::cout << "Center Theta:               " << photocathode_theta()        << std::endl;
      std::cout << "Surface Area:               " << photocathode_area           << std::endl; }
    if(verbosity > 1){
      /* placeholder */ }

    std::cout << std::endl;
  }

};


#endif
