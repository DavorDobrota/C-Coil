
#include "Polynomial.h"

void Polynomial::genDerivativeMatrix()
{
    if (!initialised)
    {
        initialised = true;
        derivativeMatrix.resize(size);

        for (int i = 0; i < size; ++i)
        {
            derivativeMatrix[i].resize(size);

            for (int j = 0; j < size; ++j)
            {
                if (j == i + 1)
                {
                    derivativeMatrix[i][j] = j;
                } else
                {
                    derivativeMatrix[i][j] = 0.0;
                }
            }
        }
    }
}

void Polynomial::setCoefficientsTo0()
{
    coefficients.resize(size);
    std::fill(coefficients.begin(), coefficients.end(), 0);
}

Polynomial::Polynomial() {

    genDerivativeMatrix();
    setCoefficientsTo0();
}

Polynomial::Polynomial(std::vector<double> &coefficientList)
{

    genDerivativeMatrix();
    setCoefficientsTo0();

    for (int i = 0; i < coefficientList.size(); ++i)
    {
        coefficients[i] = coefficientList[i];
    }
}

int Polynomial::getLeadingCoefficient()
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

std::vector<double> Polynomial::getCoefficients()
{
    std::vector<double> vector;
    int leading_coefficient = getLeadingCoefficient();

    for (int i = 0; i <= leading_coefficient; ++i)
    {
        vector.push_back(coefficients[i]);
    }
    return vector;
}

void Polynomial::printPolynomial()
{
    printf("[");
    int leading_coefficient = getLeadingCoefficient();

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

void Polynomial::printForGrapher()
{
    int leading_coefficient = getLeadingCoefficient();

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
    int leading_coefficient = getLeadingCoefficient();

    for (int i = 0; i <= leading_coefficient; ++i){
        coefficients[i] *= multiplier;
    }
}

Polynomial Polynomial::takeDerivative()
{
    Polynomial polynomial = Polynomial();

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            polynomial.coefficients[i] += derivativeMatrix[i][j] * coefficients[j];
        }
    }

    return polynomial;
}

double Polynomial::getValueAt(double x)
{
    double value = 0.0;
    int leading_coefficient = getLeadingCoefficient();

    for (int i = 0; i <= leading_coefficient; ++i)
    {
        value += coefficients[i] * pow(x, i);
    }
    return value;
}

void Polynomial::multiplyByXtoN(int n)
{
    int leading_coefficient = getLeadingCoefficient();

    if (leading_coefficient + n < size && n > 0)
    {
        for (int i = leading_coefficient; i >= 0; --i)
        {
            coefficients[i + n] = coefficients[i];
            coefficients[i] = 0.0;
        }
    }
}

Polynomial Polynomial::addPolynomials(Polynomial &pol1, Polynomial &pol2)
{
    int leading_coefficient_1 = pol1.getLeadingCoefficient();
    int leading_coefficient_2 = pol2.getLeadingCoefficient();

    Polynomial polynomial = Polynomial();

    for (int i = 0; i <= std::max(leading_coefficient_1, leading_coefficient_2); ++i)
    {
        polynomial.coefficients[i] = pol1.coefficients[i] + pol2.coefficients[i];
    }
    return polynomial;
}

Polynomial Polynomial::multiplyPolynomials(Polynomial &pol1, Polynomial &pol2)
{
    int leading_coefficient_1 = pol1.getLeadingCoefficient();
    int leading_coefficient_2 = pol2.getLeadingCoefficient();

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
    int leading_index = inputPolynomial.getLeadingCoefficient();
    inputPolynomial.multiplyByConst(1.0 / inputPolynomial.coefficients[leading_index]);

    Polynomial inputDerivative = inputPolynomial.takeDerivative();
    double currentX = firstGuess;

    for (int i = 0; i < steps; ++i)
    {
    //    printf("%.20f\n", currentX);
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

std::vector<double> Polynomial::getPolynomialRealZeros(double lowerBound, double upperBound)
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
            if (std::fabs(temp_zero - polynomial_zero) < tolerance || _isnan(temp_zero))
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

std::vector<Polynomial> Polynomial::getLegendreSequence(int maxN)
{
    std::vector<Polynomial> legendrePolynomials;

    std::vector<double> vec1;
    vec1.push_back(1);

    std::vector<double> vec2;
    vec2.push_back(0);
    vec2.push_back(1);

    legendrePolynomials.emplace_back(Polynomial(vec1));
    legendrePolynomials.emplace_back(Polynomial(vec2));

    for (int i = 2; i < maxN; ++i)
    {
        legendrePolynomials.emplace_back(
                Polynomial::genLegendrePolynomialN(i, legendrePolynomials[i-1], legendrePolynomials[i-2]));
    }
    return legendrePolynomials;
}

std::vector<double> Polynomial::getLegendreZerosForN(int n, std::vector<Polynomial> &legendreSequence)
{
    return legendreSequence[n].getPolynomialRealZeros(-1.0, 1.0);
}

std::vector<double> Polynomial::getLegendreWeightsForN(int n, std::vector<Polynomial> &legendreSequence)
{
    std::vector<double> legendre_zeros = getLegendreZerosForN(n, legendreSequence);
    std::vector<double> legendre_weight;

    for (double zero : legendre_zeros)
        {
            legendre_weight.push_back(2.0 * (1 - zero * zero) / pow((n * legendreSequence[n-1].getValueAt(zero)), 2));
        }
    return legendre_weight;
}