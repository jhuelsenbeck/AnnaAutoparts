#include <cmath>
#include <iostream>
#include "ParameterBaseFrequencies.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"

#define MIN_FREQ 0.001



ParameterBaseFrequencies::ParameterBaseFrequencies(Model* m, UserSettings* s) : Parameter(m, s, "Base Frequencies") {

    updateModifiesEigens = true;
    freqs[0].resize(4);
    freqs[1].resize(4);
    std::vector<double> alpha(4, 1.0);
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    bool err = false;
    do
        {
        err = Probability::Dirichlet::rv(&rng, alpha, freqs[0]);
        } while (anyLessThanMin(freqs[0], MIN_FREQ) == true || err == true);
    freqs[1] = freqs[0];
}

void ParameterBaseFrequencies::accept(void) {

    freqs[1] = freqs[0];
}

std::string ParameterBaseFrequencies::getHeader(void) {

    return "";
}

std::vector<double> ParameterBaseFrequencies::getValue(void) {

    std::vector<double> vals;
    for (int i=0; i<freqs[0].size(); i++)
        vals.push_back(freqs[0][i]);
    return vals;
}

std::string ParameterBaseFrequencies::getValuesAsString(int precision) {

    std::string str = "(";
    for (int i=0; i<freqs[0].size(); i++)
        {
        str += std::to_string(freqs[0][i]);
        if (i+1 != freqs[0].size())
            str += ", ";
        }
    str += ")";
    return str;
}

double ParameterBaseFrequencies::lnProbability(void) {

    return log(6.0);
}

Parameter* ParameterBaseFrequencies::newRandomInstance(void) {

    return new ParameterBaseFrequencies(modelPtr, userSettings);
}

void ParameterBaseFrequencies::print(void) {

    std::cout << getValuesAsString(3) << std::endl;
}

void ParameterBaseFrequencies::reject(void) {

    freqs[0] = freqs[1];
}

std::string ParameterBaseFrequencies::type(void) {

    return "base frequencies";
}

double ParameterBaseFrequencies::update(void) {

    UpdateInfo::updateInfo().attempt("Base Frequencies");

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double alpha0 = userSettings->getTuningBaseFrequencies();
    
    std::vector<double> oldFreqs(4);
    std::vector<double> alphaForward(4);
    for (int i=0; i<4; i++)
        {
        oldFreqs[i] = freqs[0][i];
        alphaForward[i] = freqs[0][i] * alpha0;
        }
        
    std::vector<double> newFreqs(4);
    bool err = false;
    do
        {
        err = Probability::Dirichlet::rv(&rng, alphaForward, newFreqs);
        } while (anyLessThanMin(newFreqs, MIN_FREQ) == true || err == true);
    
    std::vector<double> alphaReverse(4);
    for (int i=0; i<4; i++)
        {
        alphaReverse[i] = newFreqs[i] * alpha0;
        freqs[0][i] = newFreqs[i];
        }
        
    double lnHastingsRatio = Probability::Dirichlet::lnPdf(alphaReverse, oldFreqs) - Probability::Dirichlet::lnPdf(alphaForward, newFreqs);
    return lnHastingsRatio;
}
