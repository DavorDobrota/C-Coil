
#ifndef GENERAL_COIL_PROGRAM_POLYNOMIAL_H
#define GENERAL_COIL_PROGRAM_POLYNOMIAL_H

#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace
{
    std::vector<std::vector<double>> derivativeMatrix;
    bool initialised;
}

class Polynomial
        {
    private:
        static const int size = 100;
        static const int steps = 2000;
        constexpr static const double tolerance = 1e-9;

        std::vector<double> coefficients;

        static void genDerivativeMatrix();

        void setCoefficientsTo0();

    public:
        Polynomial();
        explicit Polynomial(std::vector<double> &coefficientList);

        int getLeadingCoefficient();
        std::vector<double> getCoefficients ();
        double getValueAt(double x);
        void printPolynomial();
        void printForGrapher();

        void multiplyByConst (double multiplier);
        void multiplyByXtoN (int n);
        Polynomial takeDerivative();

        static Polynomial addPolynomials (Polynomial &pol1, Polynomial &pol2);
        static Polynomial multiplyPolynomials (Polynomial &pol1, Polynomial &pol2);
        double findNewtonZero (double firstGuess);
        std::vector<double> getPolynomialRealZeros(double lowerBound, double upperBound);

        static Polynomial genLegendrePolynomialN (int n, Polynomial legendreN_minus_1, Polynomial legendreN_minus_2);
        static std::vector<Polynomial> getLegendreSequence (int maxN);
};

#endif //GENERAL_COIL_PROGRAM_POLYNOMIAL_H
