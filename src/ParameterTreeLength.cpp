#include <cmath>
#include <iostream>
#include "ParameterTreeLength.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"

#define MIN_LENGTH 0.001
#define MAX_LENGTH 50.0
#undef USE_GAMMA_DIRICHLET



ParameterTreeLength::ParameterTreeLength(const ParameterTreeLength& parm) : Parameter(parm) {

    std::cout << "ParameterTreeLength copy constructor" << std::endl;

    numBranches = parm.numBranches;
    alphaT = parm.alphaT;
    betaT = parm.betaT;
    lambda = parm.lambda;
    length[0] = parm.length[0];
    length[1] = parm.length[1];
}

ParameterTreeLength::ParameterTreeLength(Model* m, UserSettings* s, double lam, int nb) : Parameter(m, s, "Tree Length") {

    updateModifiesEigens = false;
    numBranches = nb;
    alphaT = 1.0;
    betaT = 1.0;
    lambda = lam;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    do
        {
        length[0] = Probability::Gamma::rv(&rng, alphaT, betaT);
        } while (length[0] < MIN_LENGTH || length[0] > MAX_LENGTH);
    length[1] = length[0];
}

ParameterTreeLength::ParameterTreeLength(Model* m, UserSettings* s, double a, double b) : Parameter(m, s, "Tree Length") {

    updateModifiesEigens = false;
    numBranches = 1.0;
    alphaT = a;
    betaT = b;
    lambda = 1.0;
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    do
        {
        length[0] = Probability::Gamma::rv(&rng, alphaT, betaT);
        } while (length[0] < MIN_LENGTH || length[0] > MAX_LENGTH);
    length[1] = length[0];
}

void ParameterTreeLength::accept(void) {

    length[1] = length[0];
}

std::string ParameterTreeLength::getHeader(void) {

    return "";
}

std::vector<double> ParameterTreeLength::getValue(void) {

    std::vector<double> vals;
    vals.push_back(length[0]);
    return vals;
}

std::string ParameterTreeLength::getValuesAsString(int precision) {

    std::string str = "";
    str += std::to_string(length[0]);
    return str;
}

double ParameterTreeLength::lnProbability(void) {

#   if defined(USE_GAMMA_DIRICHLET)
    double c = Probability::Gamma::cdf(alphaT, betaT, MAX_LENGTH) - Probability::Gamma::cdf(alphaT, betaT, MIN_LENGTH);
    return (Probability::Gamma::lnPdf(alphaT, betaT, length[0]) - log(c));
#   else
    return Probability::Gamma::lnPdf(numBranches, lambda, length[0]);
#   endif
}

Parameter* ParameterTreeLength::newRandomInstance(void) {

    return new ParameterTreeLength(modelPtr, userSettings, alphaT, betaT);
}

void ParameterTreeLength::print(void) {

    std::cout << getValuesAsString(3) << std::endl;
}

void ParameterTreeLength::reject(void) {

    length[0] = length[1];
}

std::string ParameterTreeLength::type(void) {

    return "tree length";
}

double ParameterTreeLength::update(void) {

    UpdateInfo::updateInfo().attempt("Tree Length");

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double tuning = userSettings->getTuningTreeLength();
    double newVal = 0.0, randomFactor = 1.0;
    do
        {
        randomFactor = exp( (rng.uniformRv()-0.5)*tuning );
        newVal = length[0] * randomFactor;
        } while (newVal < MIN_LENGTH || newVal > MAX_LENGTH);
    length[0] = newVal;
    return log(randomFactor);
}
