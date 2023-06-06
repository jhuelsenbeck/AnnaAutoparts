#include <cmath>
#include <iostream>
#include "ParameterExchangabilityRates.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"

#define MIN_RATES 0.001



ParameterExchangabilityRates::ParameterExchangabilityRates(const ParameterExchangabilityRates& parm) : Parameter(parm) {

    std::cout << "ParameterExchangabilityRates copy constructor" << std::endl;

    rates[0] = parm.rates[0];
    rates[1] = parm.rates[1];
}

ParameterExchangabilityRates::ParameterExchangabilityRates(Model* m, UserSettings* s) : Parameter(m, s, "Exchangeability Rates") {

    updateModifiesEigens = true;
    rates[0].resize(6);
    rates[1].resize(6);
    std::vector<double> alpha(6, 1.0);
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    bool err = false;
    do
        {
        err = Probability::Dirichlet::rv(&rng, alpha, rates[0]);
        } while (anyLessThanMin(rates[0], MIN_RATES) == true || err == true);
    rates[1] = rates[0];
}

void ParameterExchangabilityRates::accept(void) {

    rates[1] = rates[0];
}

std::string ParameterExchangabilityRates::getHeader(void) {

    return "";
}

std::vector<double> ParameterExchangabilityRates::getValue(void) {

    std::vector<double> vals;
    for (int i=0; i<rates[0].size(); i++)
        vals.push_back(rates[0][i]);
    return vals;
}

std::string ParameterExchangabilityRates::getValuesAsString(int precision) {

    std::string str = "(";
    for (int i=0; i<rates[0].size(); i++)
        {
        str += std::to_string(rates[0][i]);
        if (i+1 != rates[0].size())
            str += ", ";
        }
    str += ")";
    return str;
}

double ParameterExchangabilityRates::lnProbability(void) {

    return log(120.0);
}

Parameter* ParameterExchangabilityRates::newRandomInstance(void) {

    return new ParameterExchangabilityRates(modelPtr, userSettings);
}

void ParameterExchangabilityRates::print(void) {

    std::cout << getValuesAsString(3) << std::endl;
}

void ParameterExchangabilityRates::reject(void) {

    rates[0] = rates[1];
}

std::string ParameterExchangabilityRates::type(void) {

    return "exchangability rates";
}

double ParameterExchangabilityRates::update(void) {

    UpdateInfo::updateInfo().attempt("Exchangability Rates");

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double alpha0 = userSettings->getTuningExchangabilityRates();
    
    std::vector<double> oldRates(6);
    std::vector<double> alphaForward(6);
    for (int i=0; i<6; i++)
        {
        oldRates[i] = rates[0][i];
        alphaForward[i] = rates[0][i] * alpha0;
        }
        
    std::vector<double> newRates(6);
    bool err = false;
    do
        {
        err = Probability::Dirichlet::rv(&rng, alphaForward, newRates);
        } while (anyLessThanMin(newRates, MIN_RATES) == true || err == true);
    
    std::vector<double> alphaReverse(6);
    for (int i=0; i<6; i++)
        {
        alphaReverse[i] = newRates[i] * alpha0;
        rates[0][i] = newRates[i];
        }
        
    double lnHastingsRatio = Probability::Dirichlet::lnPdf(alphaReverse, oldRates) - Probability::Dirichlet::lnPdf(alphaForward, newRates);
    return lnHastingsRatio;
}
