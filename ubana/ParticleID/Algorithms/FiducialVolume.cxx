#include "FiducialVolume.h"

namespace fidvol{

  std::vector<double> fiducialVolume::setFiducialVolume(std::vector<double> fv, fhicl::ParameterSet const & p){

    double xl = p.get< double > ("X_LOW", 10);
    double xh = p.get< double > ("X_HIGH", 10);
    double yl = p.get< double > ("Y_LOW", 10);
    double yh = p.get< double > ("Y_HIGH", 10);
    double zl = p.get< double > ("Z_LOW", 10);
    double zh = p.get< double > ("Z_HIGH", 10);

    fv.push_back(xl);
    fv.push_back(xh);
    fv.push_back(yl);
    fv.push_back(yh);
    fv.push_back(zl);
    fv.push_back(zh);

    return fv;

  }

  void fiducialVolume::printFiducialVolume(std::vector<double> fv){

    std::cout << "----- PRINTING FIDUCIAL VOLUME INFORMATION -----" << std::endl;
    std::cout << "X_LOW: " << fv.at(0) << std::endl;
    std::cout << "X_HIGH: " << fv.at(1) << std::endl;
    std::cout << "Y_LOW: " << fv.at(2) << std::endl;
    std::cout << "Y_HIGH: " << fv.at(3) << std::endl;
    std::cout << "Z_LOW: " << fv.at(4) << std::endl;
    std::cout << "Z_HIGH: " << fv.at(5) << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

  }

  bool fiducialVolume::isInFiducialVolume(TVector3 vec, std::vector<double> fv){

    double x = vec.X();
    double y = vec.Y();
    double z = vec.Z();

    std::vector<double> tpcd = fiducialVolume::getTpcDimensions();

    if (x > (tpcd.at(0)+fv.at(0)) && x < (tpcd.at(1)-fv.at(1)) && 
        y > (tpcd.at(2)+fv.at(2)) && y < (tpcd.at(3)-fv.at(3)) && 
        z > (tpcd.at(4)+fv.at(4)) && z < (tpcd.at(5)-fv.at(5))) 
      return true;

    else return false;

  }

  std::vector<double> fiducialVolume::getTpcDimensions(){

    std::vector<double> tpcd;
    double tpc_xl = 0.0;
    double tpc_xh = 256.0;
    double tpc_yl = -116.5;
    double tpc_yh = 116.5;
    double tpc_zl = 0.0;
    double tpc_zh = 1036.0;

    tpcd.push_back(tpc_xl);
    tpcd.push_back(tpc_xh);
    tpcd.push_back(tpc_yl);
    tpcd.push_back(tpc_yh);
    tpcd.push_back(tpc_zl);
    tpcd.push_back(tpc_zh);

    return tpcd;
  }

}
