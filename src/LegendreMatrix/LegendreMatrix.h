#ifndef GENERAL_COIL_PROGRAM_LEGENDREMATRIX_H
#define GENERAL_COIL_PROGRAM_LEGENDREMATRIX_H

/// @brief Contains matrices with precomputed Gauss-Legendre quadrature weights and positions up to maxLegendreOrder.
/// Row defines the quadrature order n, and column the appropriate index i (up to n).
namespace Legendre
{
    const int maxLegendreOrder = 100;

    const extern double positionMatrix[maxLegendreOrder][maxLegendreOrder];
    const extern double weightsMatrix[maxLegendreOrder][maxLegendreOrder];
}

#endif //GENERAL_COIL_PROGRAM_LEGENDREMATRIX_H
