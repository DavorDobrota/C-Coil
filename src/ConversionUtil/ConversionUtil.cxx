#include <cmath>
#include <vector>

#include "ConversionUtil.h"

void convertPolarToCylindrical(double polarR, double polarTheta, double polarPhi,
                                     double &cylindricalZ, double &cylindricalR, double &cylindricalPhi)
{
    cylindricalZ = polarR * cos(polarTheta);
    cylindricalR = polarR * sin(polarTheta);
    cylindricalPhi = cylindricalPhi;
}

void convertAllPolarToCylindrical(const std::vector<double> &polarRArr,
                                        const std::vector<double> &polarThetaArr,
                                        const std::vector<double> &polarPhiArr,
                                        std::vector<double> &cylindricalZArr,
                                        std::vector<double> &cylindricalRArr,
                                        std::vector<double> &cylindricalPhiArr)
{
    if (polarRArr.size() == polarThetaArr.size() == polarPhiArr.size())
    {
        cylindricalZArr.resize(0);
        cylindricalRArr.resize(0);
        cylindricalPhiArr.resize(0);

        for (int i = 0; i < polarRArr.size(); ++i)
        {
            cylindricalZArr.push_back(polarRArr[i] * cos(polarThetaArr[i]));
            cylindricalRArr.push_back(polarRArr[i] * sin(polarThetaArr[i]));
            cylindricalPhiArr.push_back(polarPhiArr[i]);
        }
    }
}

void convertCylindricalToPolar(double cylindricalZ, double cylindricalR, double cylindricalPhi,
                                     double &polarR, double &polarTheta, double &polarPhi)
{
    polarR = std::sqrt(cylindricalZ * cylindricalZ + cylindricalR * cylindricalR);
    polarTheta = atan2(cylindricalR, cylindricalZ);
    polarPhi = cylindricalPhi;
}

void convertAllCylindricalToPolar(const std::vector<double> &cylindricalZArr,
                                        const std::vector<double> &cylindricalRArr,
                                        const std::vector<double> &cylindricalPhiArr,
                                        std::vector<double> &polarRArr,
                                        std::vector<double> &polarThetaArr,
                                        std::vector<double> &polarPhiArr)
{
    if (cylindricalZArr.size() == cylindricalRArr.size() && cylindricalRArr.size() == cylindricalPhiArr.size())
    {
        polarRArr.resize(0);
        polarThetaArr.resize(0);
        polarPhiArr.resize(0);

        for (int i = 0; i < cylindricalZArr.size(); ++i)
        {
            polarRArr.push_back(std::sqrt(cylindricalZArr[i] * cylindricalZArr[i] + cylindricalRArr[i] * cylindricalRArr[i]));
            polarThetaArr.push_back(atan2(cylindricalRArr[i], cylindricalZArr[i]));
            polarPhiArr.push_back(cylindricalPhiArr[i]);
        }
    }
}
