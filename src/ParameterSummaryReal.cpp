#include <iomanip>
#include <iostream>
#include <algorithm>
#include "Model.hpp"
#include "ParameterSummaryReal.hpp"



std::pair<double,double> ParameterSummaryReal::calculateCredibleInterval(double cv) {

    std::sort(values.begin(), values.end());
    size_t size = values.size();
    size_t lp = size / 40;
    size_t up = size - lp - 1;
    std::pair<double,double> ci = std::make_pair(values[lp], values[up]);
    return ci;
}

double ParameterSummaryReal::calculateMean(void) {

    double a = 0.0;
    for (int i=0; i<values.size(); i++)
        {
        double x = values[i];
        if (i == 0)
            {
            a = x;
            }
        else
            {
            double aOld = a;
            a = aOld + (x - aOld) / (i+1);
            }
        }

    double mean = a;
        
    return mean;
}

double ParameterSummaryReal::calculateVariance(void) {

    double a = 0.0, s = 0.0;
    for (int i=0; i<values.size(); i++)
        {
        double x = values[i];
        if (i == 0)
            {
            a = x;
            s = 0.0;
            }
        else
            {
            double aOld = a;
            a = aOld + (x - aOld) / (i+1);
            s = s + (x - aOld) * (x - aOld);
            }
        }

    double variance = 0.0;
    if (values.size() > 1)
        variance = s / (values.size() - 1);
        
    return variance;
}

double ParameterSummaryReal::getTrueValue(Model* m) {

    // find the restaurant for the parameter
    std::vector<std::string> nameTokens;
    std::string str = "";
    for (int i=0; i<name.length(); i++)
        {
        char c = name[i];
        if (c != ',' && c != '(' && c != ')')
            str += std::string(1,c);
        else
            {
            if (str != "")
                nameTokens.push_back(str);
            str = "";
            }
        }
    
    double trueValue = 0.0;
    if ( !(nameTokens.size() == 0 || nameTokens[1] == "Alpha") )
        {
        int subset = stoi(nameTokens[1]);
        if (nameTokens[0] == "L")
            {
            // tree length restaurant
            trueValue = m->getTreeLength(subset-1);
            }
        else if (nameTokens[0] == "Pi")
            {
            // base frequency restaurant
            std::vector<double>& f = m->getBaseFrequencies(subset-1);
            if (nameTokens[2] == "A")
                trueValue = f[0];
            else if (nameTokens[2] == "C")
                trueValue = f[1];
            else if (nameTokens[2] == "G")
                trueValue = f[2];
            else
                trueValue = f[3];
            }
        else if (nameTokens[0] == "R")
            {
            // exchangability rates restaurant
            std::vector<double>& r = m->getExchangeabilityRates(subset-1);
            if (nameTokens[2] == "AC")
                trueValue = r[0];
            else if (nameTokens[2] == "AG")
                trueValue = r[1];
            else if (nameTokens[2] == "AT")
                trueValue = r[2];
            else if (nameTokens[2] == "CG")
                trueValue = r[3];
            else if (nameTokens[2] == "CT")
                trueValue = r[4];
            else
                trueValue = r[5];
            }
        else if (nameTokens[0] == "Alpha")
            {
            // gamma shape restaurant
            trueValue = m->getGammaShape(subset-1);
            }
        }
    else if (nameTokens.size() != 0 && nameTokens[1] == "Alpha")
        {
        if (nameTokens[0] == "L")
            trueValue = m->getAlphaTreeLength();
        else if (nameTokens[0] == "Pi")
            trueValue = m->getAlphaFrequencies();
        else if (nameTokens[0] == "R")
            trueValue = m->getAlphaRates();
        else if (nameTokens[0] == "Alpha")
            trueValue = m->getAlphaShape();
        }
        
    return trueValue;
}

bool ParameterSummaryReal::inCi(Model* m) {

    double trueValue = getTrueValue(m);
    std::pair<double,double> ci = calculateCredibleInterval(0.95);
    bool inCi = true;
    if (trueValue < ci.first || trueValue > ci.second)
        inCi = false;
    return inCi;
}

double ParameterSummaryReal::mse(Model* m) {

    double trueValue = getTrueValue(m);

    double a = 0.0;
    for (int i=0; i<values.size(); i++)
        {
        double x = values[i];
        a += (trueValue - x) * (trueValue - x);
        }
    a /= values.size();

    return a;
}

void ParameterSummaryReal::print(void) {
    
    std::cout << "   * " << name << ": ";
    double mean = calculateMean();
    double variance = calculateVariance();
    std::pair<double,double> ci = calculateCredibleInterval(0.95);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << mean << " " << variance << " (" << ci.first << "," << ci.second << ")" << std::endl;
}

void ParameterSummaryReal::print(Model* m) {
        
    double trueValue = getTrueValue(m);

    std::cout << "   * " << name << ": ";
    double mean = calculateMean();
    double variance = calculateVariance();
    std::pair<double,double> ci = calculateCredibleInterval(0.95);
    bool inCi = true;
    if (trueValue < ci.first || trueValue > ci.second)
        inCi = false;
    double mse = (trueValue - mean) * (trueValue - mean);
    std::cout << std::fixed << std::setprecision(4);
    std::cout << mean << " " << variance << " (" << ci.first << "," << ci.second << ") (" << trueValue << " " << mse << ((inCi == true) ? " YES" : " NO") << ")" << std::endl;
}

std::string ParameterSummaryReal::summarize(void) {
    
    double mean = calculateMean();
    double variance = calculateVariance();
    std::pair<double,double> ci = calculateCredibleInterval(0.95);
    
    std::string str = "";
    str += std::to_string(mean);
    str += " ";
    str += std::to_string(variance);
    str += " (";
    str += std::to_string(ci.first);
    str += ", ";
    str += std::to_string(ci.second);
    str += ")";

    return str;
}
