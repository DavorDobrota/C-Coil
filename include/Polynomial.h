
#ifndef GENERAL_COIL_PROGRAM_POLYNOMIAL_H
#define GENERAL_COIL_PROGRAM_POLYNOMIAL_H

#include <vector>
#include <cstdio>
#include <cmath>

namespace
{
    std::vector<std::vector<double>> derivativeMatrix;
    bool initialised;
}

class Polynomial
        {
    private:
        static const int size = 100;
        static const int steps = 10;

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

        void multiplyByConst (double multiplier);
        void multiplyByXtoN (int n);
        Polynomial takeDerivative();

        static Polynomial addPolynomials (Polynomial pol1, Polynomial pol2);
        static Polynomial multiplyPolynomials (Polynomial pol1, Polynomial pol2);
        static double findNewtonZero (double firstGuess, Polynomial &inputPolynomial);
        static Polynomial genLegendrePolynomialN (int n, Polynomial legendreN_minus_1, Polynomial legendreN_minus_2);

};


#endif //GENERAL_COIL_PROGRAM_POLYNOMIAL_H
