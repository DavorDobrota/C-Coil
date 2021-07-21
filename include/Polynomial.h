#ifndef GENERAL_COIL_PROGRAM_POLYNOMIAL_H
#define GENERAL_COIL_PROGRAM_POLYNOMIAL_H

#include <vector>

class Polynomial
{
private:
    constexpr static int size = 25;
    constexpr static int steps = 50;
    constexpr static double tolerance = 1e-10;

    std::vector<double> coefficients;

    static void genDerivativeMatrix();

    void setCoefficientsToZero();

public:
    Polynomial();
    explicit Polynomial(const std::vector<double> &coefficientList);

    int getLeadingCoefficientIndex() const;
    std::vector<double> getCoefficients() const;
    double getValueAt(double x) const;
    void printPolynomial() const;
    void printForGrapher() const;

    void multiplyByConst(double multiplier);
    void multiplyByXtoN(int n);
    Polynomial takeDerivative() const;

    static Polynomial addPolynomials(const Polynomial &pol1, const Polynomial &pol2);
    static Polynomial multiplyPolynomials(const Polynomial &pol1, const Polynomial &pol2);
    static double findNewtonZero(double firstGuess, Polynomial inputPolynomial);
    static double findHouseholderZero(double firstGuess, Polynomial inputPolynomial);
    std::vector<double> getPolynomialRealZeros(double lowerBound, double upperBound) const;

    static Polynomial genLegendrePolynomialN(int n, Polynomial legendreN_minus_1, Polynomial legendreN_minus_2);
    static std::vector<Polynomial> getLegendreSequenceUpToN(int maxN);
    static void genInternalLegendreSequence();
    static void getLegendreParametersForN(int n, std::vector<double> &zeros, std::vector<double> &weights);
};

#endif //GENERAL_COIL_PROGRAM_POLYNOMIAL_H
