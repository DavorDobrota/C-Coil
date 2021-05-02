
#ifndef GENERAL_COIL_PROGRAM_POLYNOMIAL_H
#define GENERAL_COIL_PROGRAM_POLYNOMIAL_H

#include <vector>


class Polynomial
{
private:
    constexpr static int size = 100;
    constexpr static int steps = 50;
    constexpr static double tolerance = 1e-9;

    std::vector<double> coefficients;

    static void genDerivativeMatrix();

    void setCoefficientsToZero();

public:
    Polynomial();
    explicit Polynomial(std::vector<double> &coefficientList);

    int getLeadingCoefficientIndex();
    std::vector<double> getCoefficients();
    double getValueAt(double x);
    void printPolynomial();
    void printForGrapher();

    void multiplyByConst(double multiplier);
    void multiplyByXtoN(int n);
    Polynomial takeDerivative();

    static Polynomial addPolynomials(Polynomial &pol1, Polynomial &pol2);
    static Polynomial multiplyPolynomials(Polynomial &pol1, Polynomial &pol2);
    static double findNewtonZero(double firstGuess, Polynomial inputPolynomial);
    static double findHouseholderZero(double firstGuess, Polynomial inputPolynomial);
    std::vector<double> getPolynomialRealZeros(double lowerBound, double upperBound);

    static Polynomial genLegendrePolynomialN(int n, Polynomial legendreN_minus_1, Polynomial legendreN_minus_2);
    static std::vector<Polynomial> getLegendreSequenceUpToN(int maxN);
    static void genInternalLegendreSequence();
    static void getLegendreParametersForN(int n, std::vector<double> &zeros, std::vector<double> &weights);
};

#endif //GENERAL_COIL_PROGRAM_POLYNOMIAL_H
