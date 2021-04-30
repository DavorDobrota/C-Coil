
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

Polynomial Polynomial::addPolynomials(Polynomial pol1, Polynomial pol2)
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

double Polynomial::findNewtonZero(double firstGuess, Polynomial &inputPolynomial)
{
    Polynomial inputDerivative = inputPolynomial.takeDerivative();
    double currentX = firstGuess;

    for (int i = 0; i < steps; ++i)
    {
        currentX -= inputPolynomial.getValueAt(currentX) / inputDerivative.getValueAt(currentX);
    }
    return currentX;
}

Polynomial Polynomial::genLegendrePolynomialN(int n, Polynomial legendreN_minus_1, Polynomial legendreN_minus_2)
{
    legendreN_minus_1.multiplyByXtoN(1);
    legendreN_minus_1.multiplyByConst((2 * n - 1));
    legendreN_minus_2.multiplyByConst(n - 1);

    Polynomial polynomial = addPolynomials(legendreN_minus_1, legendreN_minus_2);
    polynomial.multiplyByConst(1.0 / n);

    return polynomial;
}