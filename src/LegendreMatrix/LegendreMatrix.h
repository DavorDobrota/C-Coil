#ifndef GENERAL_COIL_PROGRAM_LEGENDREMATRIX_H
#define GENERAL_COIL_PROGRAM_LEGENDREMATRIX_H

namespace Legendre
{
    const int maxLegendreOrder = 64;

    const extern double positionMatrix[maxLegendreOrder][maxLegendreOrder];
    const extern double weightsMatrix[maxLegendreOrder][maxLegendreOrder];
}

#endif //GENERAL_COIL_PROGRAM_LEGENDREMATRIX_H
