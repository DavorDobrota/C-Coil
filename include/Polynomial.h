
#ifndef GENERAL_COIL_PROGRAM_POLYNOMIAL_H
#define GENERAL_COIL_PROGRAM_POLYNOMIAL_H

#include <vector>
#include <string>
#include <cmath>

class Polynomial {

private:
    static const int size = 100;
    static const int steps = 10;

    static std::vector<std::vector<double>> derivativeMatrix;
    static bool initialised;

    std::vector<double> coefficients;

    static void genDerivativeMatrix();

    void setCoefficientsTo0();

public:
    Polynomial();
    Polynomial(int degree, std::vector<double> &coefficientList);

    int getLeadingCoefficient();
    std::vector<double> getCoefficients ();
    double getValueAt(double x);
    std::string toString();

    static Polynomial derivativeOfPolynomial(Polynomial &inputPolynomial);
    static double findNewtonZero (double firstGuess, Polynomial &inputPolynomial);
    static Polynomial genLegendrePolynomialN (int n, Polynomial &legendreN_minus_1, Polynomial &legendreN_minus_2);

};


#endif //GENERAL_COIL_PROGRAM_POLYNOMIAL_H
