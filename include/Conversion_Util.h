
#ifndef GENERAL_COIL_PROGRAM_CONVERSION_UTIL_H
#define GENERAL_COIL_PROGRAM_CONVERSION_UTIL_H

#include <vector>

void convertPolarToCylindrical(double polarR, double polarTheta, double polarPhi,
                               double &cylindricalZ, double &cylindricalR, double &cylindricalPhi);

void convertAllPolarToCylindrical(const std::vector<double> &polarRArr,
                                  const std::vector<double> &polarThetaArr,
                                  const std::vector<double> &polarPhiArr,
                                  std::vector<double> &cylindricalZArr,
                                  std::vector<double> &cylindricalRArr,
                                  std::vector<double> &cylindricalPhiArr);

void convertCylindricalToPolar(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                               double &polarR, double &polarTheta, double &polarPhi);

void convertAllCylindricalToPolar(const std::vector<double> &cylindricalZArr,
                                  const std::vector<double> &cylindricalRArr,
                                  const std::vector<double> &cylindricalPhiArr,
                                  std::vector<double> &polarRArr,
                                  std::vector<double> &polarThetaArr,
                                  std::vector<double> &polarPhiArr);

// TODO - add Cartesian transformations

#endif //GENERAL_COIL_PROGRAM_CONVERSION_UTIL_H
