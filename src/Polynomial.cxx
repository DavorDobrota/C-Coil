
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

Polynomial::Polynomial(int degree, std::vector<double> &coefficientList)
{

    genDerivativeMatrix();
    setCoefficientsTo0();

    for (int i = 0; i <= degree; ++i)
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

std::string Polynomial::toString()
{
    std::string output = "[";
    int leading_coefficient = getLeadingCoefficient();

    for (int i = 0; i <= leading_coefficient; ++i)
    {
        output += coefficients[i];
        if (i != leading_coefficient)
        {
            output += ", ";
        }
    }
    return output;
}

Polynomial Polynomial::derivativeOfPolynomial(Polynomial &inputPolynomial) {

    Polynomial polynomial = Polynomial();

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            polynomial.coefficients[i] += derivativeMatrix[i][j] * inputPolynomial.coefficients[j];
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

double Polynomial::findNewtonZero(double firstGuess, Polynomial &inputPolynomial)
{
    Polynomial inputDerivative = derivativeOfPolynomial(inputPolynomial);
    double currentX = firstGuess;

    for (int i = 0; i < steps; ++i)
    {
        currentX -= inputPolynomial.getValueAt(currentX) / inputDerivative.getValueAt(currentX);
    }
    return currentX;
}