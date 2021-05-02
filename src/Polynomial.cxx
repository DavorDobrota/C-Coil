#include <cstdio>
#include <cmath>

#include <algorithm>
#include <vector>

#include "Polynomial.h"

namespace
{
    std::vector<std::vector<double>> g_derivativeMatrix;
    bool g_initialised;
    std::vector<Polynomial> g_legendreSequence;
    bool g_legendreInitialised;
}

void Polynomial::genDerivativeMatrix()
{
    if (!g_initialised)
    {
        g_initialised = true;
        g_derivativeMatrix.resize(size);

        for (int i = 0; i < size; ++i)
        {
            g_derivativeMatrix[i].resize(size);

            for (int j = 0; j < size; ++j)
            {
                if (j == i + 1)
                {
                    g_derivativeMatrix[i][j] = j;
                } else
                {
                    g_derivativeMatrix[i][j] = 0.0;
                }
            }
        }
    }
}

void Polynomial::setCoefficientsToZero()
{
    coefficients.resize(size);
    std::fill(coefficients.begin(), coefficients.end(), 0);
}

Polynomial::Polynomial()
{
    genDerivativeMatrix();
    setCoefficientsToZero();
}

Polynomial::Polynomial(const std::vector<double> &coefficientList)
{

    genDerivativeMatrix();
    setCoefficientsToZero();

    for (int i = 0; i < coefficientList.size(); ++i)
    {
        coefficients[i] = coefficientList[i];
    }
}

int Polynomial::getLeadingCoefficientIndex() const
{
    int leading_coefficient;

    for (int i = size - 1; i >= 0; --i)
    {
        if (coefficients[i] != 0)
        {
            leading_coefficient = i;
            break;
        }
    }
    return leading_coefficient;
}

std::vector<double> Polynomial::getCoefficients() const
{
    int leading_coefficient = getLeadingCoefficientIndex();

    return std::vector<double>(coefficients.begin(), coefficients.begin() + leading_coefficient + 1);
}

void Polynomial::printPolynomial() const
{
    printf("[");
    int leading_coefficient = getLeadingCoefficientIndex();

    for (int i = 0; i <= leading_coefficient; ++i)
    {
        printf("%.4f", coefficients[i]);
        if (i != leading_coefficient)
        {
            printf(", ");
        }
    }
    printf("]\n");
}

void Polynomial::printForGrapher() const
{
    int leading_coefficient = getLeadingCoefficientIndex();

    for (int i = 0; i <= leading_coefficient; ++i)
    {
        if (std::abs(coefficients[i]) > tolerance)
        {
            printf("%fx^%d", coefficients[i], i);
            if (i != leading_coefficient)
            {
                printf(" + ");
            }
        }
    }
    printf("\n");
}

void Polynomial::multiplyByConst(double multiplier)
{
    int leading_coefficient = getLeadingCoefficientIndex();

    for (int i = 0; i <= leading_coefficient; ++i){
        coefficients[i] *= multiplier;
    }
}

Polynomial Polynomial::takeDerivative() const
{
    Polynomial polynomial = Polynomial();

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            polynomial.coefficients[i] += g_derivativeMatrix[i][j] * coefficients[j];
        }
    }

    return polynomial;
}

double Polynomial::getValueAt(double x) const
{
    double value = 0.0;
    int leading_coefficient = getLeadingCoefficientIndex();

    for (int i = 0; i <= leading_coefficient; ++i)
    {
        value += coefficients[i] * pow(x, i);
    }
    return value;
}

void Polynomial::multiplyByXtoN(int n)
{
    int leading_coefficient = getLeadingCoefficientIndex();

    if (leading_coefficient + n < size && n > 0)
    {
        for (int i = leading_coefficient; i >= 0; --i)
        {
            coefficients[i + n] = coefficients[i];
            coefficients[i] = 0.0;
        }
    }
}

Polynomial Polynomial::addPolynomials(const Polynomial &pol1, const Polynomial &pol2)
{
    int leading_coefficient_1 = pol1.getLeadingCoefficientIndex();
    int leading_coefficient_2 = pol2.getLeadingCoefficientIndex();

    Polynomial polynomial = Polynomial();

    for (int i = 0; i <= std::max(leading_coefficient_1, leading_coefficient_2); ++i)
    {
        polynomial.coefficients[i] = pol1.coefficients[i] + pol2.coefficients[i];
    }
    return polynomial;
}

Polynomial Polynomial::multiplyPolynomials(const Polynomial &pol1, const Polynomial &pol2)
{
    int leading_coefficient_1 = pol1.getLeadingCoefficientIndex();
    int leading_coefficient_2 = pol2.getLeadingCoefficientIndex();

    Polynomial polynomial = Polynomial();

    if (leading_coefficient_1 + leading_coefficient_2 < size)
    {
        for (int i = 0; i <= leading_coefficient_1; ++i)
        {
            for (int j = 0; j <= leading_coefficient_2; ++j)
            {
                polynomial.coefficients[i+j] += pol1.coefficients[i] * pol2.coefficients[j];
            }
        }
    }
    return polynomial;
}

double Polynomial::findNewtonZero(double firstGuess, Polynomial inputPolynomial)
{
    int leading_index = inputPolynomial.getLeadingCoefficientIndex();
    inputPolynomial.multiplyByConst(1.0 / inputPolynomial.coefficients[leading_index]);

    Polynomial inputDerivative = inputPolynomial.takeDerivative();
    double currentX = firstGuess;

    for (int i = 0; i < steps; ++i)
    {
        currentX -= inputPolynomial.getValueAt(currentX) / inputDerivative.getValueAt(currentX);

    }
    return currentX;
}

double Polynomial::findHouseholderZero(double firstGuess, Polynomial inputPolynomial)
{
    Polynomial inputFirstDerivative = inputPolynomial.takeDerivative();
    Polynomial inputSecondDerivative = inputFirstDerivative.takeDerivative();
    double currentX = firstGuess;

    double temp_h;

    for (int i = 0; i < steps; ++i)
    {
        temp_h = inputPolynomial.getValueAt(currentX) / inputFirstDerivative.getValueAt(currentX);
        currentX -= temp_h / (1 - 0.5 * temp_h *
                (inputSecondDerivative.getValueAt(currentX) / inputFirstDerivative.getValueAt(currentX)));
    }
    return currentX;
}

std::vector<double> Polynomial::getPolynomialRealZeros(double lowerBound, double upperBound) const
{
    double zeros_from_Attempts[size * 2 + 1];
    double step = (upperBound - lowerBound) / (2 * size);

    for (int i = 0; i <= 2 * size; ++i)
    {
        zeros_from_Attempts[i] = findHouseholderZero(lowerBound + i * step, *this);
    }

    std::vector<double> polynomial_zeros;
    double temp_zero;

    for (int i = 0; i <= 2 * size; ++i)
    {
        temp_zero = zeros_from_Attempts[i];
        bool duplicate = false;

        for (double polynomial_zero : polynomial_zeros)
        {
            if (std::fabs(temp_zero - polynomial_zero) < tolerance || std::isnan(temp_zero))
            {
                duplicate = true;
            }
        }
        if (!duplicate)
        {
            polynomial_zeros.push_back(temp_zero);
        }
    }
    std::sort(polynomial_zeros.begin(), polynomial_zeros.end());

    return polynomial_zeros;
}

Polynomial Polynomial::genLegendrePolynomialN(int n, Polynomial legendreN_minus_1, Polynomial legendreN_minus_2)
{
    legendreN_minus_1.multiplyByXtoN(1);
    legendreN_minus_1.multiplyByConst((2 * n - 1));
    legendreN_minus_2.multiplyByConst(- n + 1);

    Polynomial polynomial = addPolynomials(legendreN_minus_1, legendreN_minus_2);
    polynomial.multiplyByConst(1.0 / n);

    return polynomial;
}

std::vector<Polynomial> Polynomial::getLegendreSequenceUpToN(int maxN)
{
    std::vector<Polynomial> legendrePolynomials;

    std::vector<double> vec1 = {1};
    std::vector<double> vec2 = {0, 1};

    legendrePolynomials.emplace_back(Polynomial(vec1));
    legendrePolynomials.emplace_back(Polynomial(vec2));

    for (int i = 2; i < maxN; ++i)
    {
        legendrePolynomials.emplace_back(
                Polynomial::genLegendrePolynomialN(i, legendrePolynomials[i-1], legendrePolynomials[i-2]));
    }
    return legendrePolynomials;
}

void Polynomial::genInternalLegendreSequence()
{
    if (!g_legendreInitialised){
        g_legendreSequence = getLegendreSequenceUpToN(size);
        g_legendreInitialised = true;
    }
}

void Polynomial::getLegendreParametersForN(int n, std::vector<double> &zeros,
                                                  std::vector<double> &weights)
{
    genInternalLegendreSequence();
    zeros = g_legendreSequence[n].getPolynomialRealZeros(-1.0, 1.0);
    weights.resize(0);

    for (double zero : zeros)
    {
        weights.push_back(2.0 * (1 - zero * zero) / pow((n * g_legendreSequence[n-1].getValueAt(zero)), 2));
    }
}
