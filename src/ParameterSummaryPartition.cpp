#include <iomanip>
#include <iostream>
#include "Model.hpp"
#include "ParameterSummaryPartition.hpp"
#include "Partition.hpp"
#include "RandomVariable.hpp"



bool cmp(std::pair<std::string, int>& a, std::pair<std::string, int>& b) {

    return a.second > b.second;
}

ParameterSummaryPartition::ParameterSummaryPartition(void) : ParameterSummary() {

    numSamples = 0;
}

ParameterSummaryPartition::~ParameterSummaryPartition(void) {

    for (std::map<Partition*,int>::iterator it = partitionsMap.begin(); it != partitionsMap.end(); it++)
        delete it->first;
}

void ParameterSummaryPartition::addValue(std::string s) {

    std::map<std::string,int>::iterator it = partitionsStringMap.find(s);
    if (it == partitionsStringMap.end())
        {
        partitionsStringMap.insert( std::make_pair(s,1) );
        }
    else
        {
        it->second++;
        }
    numSamples++;
}

void ParameterSummaryPartition::calculateMeanPartition(void) {

    // initialize the partition map, if it does not yet exist
    if (partitionsMap.size() == 0)
        {
        for (auto& it : partitionsStringMap)
            {
            // change to vector
            std::vector<std::string> partitionTokens;
            std::string partStr = it.first;
            std::string token = "";
            for (int i=0; i<partStr.length(); i++)
                {
                char c = partStr[i];
                if (c == ',')
                    {
                    if (token != "")
                        {
                        partitionTokens.push_back(token);
                        token = "";
                        }
                    }
                else
                    token += std::string(1,c);
                }
            if (token != "")
                partitionTokens.push_back(token);
                
            std::vector<int> vec;
            for (int i=0; i<partitionTokens.size(); i++)
                {
                int x = stoi(partitionTokens[i]);
                vec.push_back(x);
                }
                
            
            Partition* part = new Partition(vec);
            partitionsMap.insert( std::make_pair(part,it.second) );
            }
        }
        
        
    //for (std::map<Partition*,int>::iterator it = partitionsMap.begin(); it != partitionsMap.end(); it++)
    //    std::cout << it->first->getRgfString() << " " << it->second << std::endl;
        
    double score = 0.0;
    Partition meanPartition(partitionsMap, score);
    std::cout << "   * Mean partition: " << meanPartition.getRgfString() << " (" << score << ")" << std::endl;
    
}

std::string ParameterSummaryPartition::getTruePartition(Model* m) {

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

    std::string truePartition = "";
    if (nameTokens[0] == "RGF")
        {
        if (nameTokens[1] == "Length")
            {
            truePartition = m->getPartitionTreeLength();
            }
        else if (nameTokens[1] == "Pi")
            {
            truePartition = m->getPartitionBaseFrequencies();
            }
        else if (nameTokens[1] == "R")
            {
            truePartition = m->getPartitionRates();
           }
        else if (nameTokens[1] == "Alpha")
            {
            truePartition = m->getPartitionGammaShape();
            }
        }
    return truePartition;
}

bool ParameterSummaryPartition::inCi(Model* m) {

    std::string truePartition = getTruePartition(m);

    std::vector<std::pair<std::string, int> > v;
    for (auto& it : partitionsStringMap)
        v.push_back(it);
    sort(v.begin(), v.end(), cmp);

    double cumulativeProb = 0.0;
    for (int i=0; i<v.size(); i++)
        {
        double prob = (double)v[i].second / numSamples;
        if (v[i].first == truePartition)
            {
            bool inCi = false;
            if (cumulativeProb + prob <= 0.95)
                inCi = true;
            else if (cumulativeProb + prob > 0.95 && cumulativeProb <= 0.95)
                {
                double u = RandomVariable::randomVariableInstance().uniformRv();
                if (u < (0.95 - cumulativeProb) / prob)
                    inCi = true;
                }
            if (inCi == true)
                return true;
            }
        cumulativeProb += prob;
        }

    return false;
}

double ParameterSummaryPartition::mse(Model* m) {

    return 0.0;
}

void ParameterSummaryPartition::pairProbabilities(void) {

    if (partitionsMap.size() == 0)
        {
        for (auto& it : partitionsStringMap)
            {
            // change to vector
            std::vector<std::string> partitionTokens;
            std::string partStr = it.first;
            std::string token = "";
            for (int i=0; i<partStr.length(); i++)
                {
                char c = partStr[i];
                if (c == ',')
                    {
                    if (token != "")
                        {
                        partitionTokens.push_back(token);
                        token = "";
                        }
                    }
                else
                    token += std::string(1,c);
                }
            if (token != "")
                partitionTokens.push_back(token);
                
            std::vector<int> vec;
            for (int i=0; i<partitionTokens.size(); i++)
                {
                int x = stoi(partitionTokens[i]);
                vec.push_back(x);
                }
                
            
            Partition* part = new Partition(vec);
            partitionsMap.insert( std::make_pair(part,it.second) );
            }
        }

    int n = 0;
    for (std::map<Partition*,int>::iterator it = partitionsMap.begin(); it != partitionsMap.end(); it++)
        n += it->second;
        
    int np = partitionsMap.begin()->first->getNumElements();
    std::cout << "   * Pair Probabilities: ";
    for (int i=0; i<np; i++)
        {
        for (int j=i+1; j<np; j++)
            {
            if ( !(i == 0 && j == i+1) )
                std::cout << "   *                     ";
            int numSame = 0;
            for (std::map<Partition*,int>::iterator it = partitionsMap.begin(); it != partitionsMap.end(); it++)
                {
                Subset* sI = it->first->findSubsetWithElement(i);
                Subset* sJ = it->first->findSubsetWithElement(j);
                if (sI == sJ)
                    numSame += it->second;
                }
            std::cout << i+1 << "-" << j+1 << " -- " << (double)numSame/n << std::endl;
            }
        }
}

void ParameterSummaryPartition::print(void) {

    std::vector<std::pair<std::string, int> > v;
    for (auto& it : partitionsStringMap)
        v.push_back(it);
    sort(v.begin(), v.end(), cmp);

    int nameSize = (int)name.length() + 2;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "   * " << name << ": ";
    double cumulativeProb = 0.0;
    double numInCi = 0.0;
    for (int i=0; i<v.size(); i++)
        {
        if (i != 0)
            {
            std::cout << "   * ";
            for (int i=0; i<nameSize; i++)
                std::cout << " ";
            }
        double prob = (double)v[i].second / numSamples;
        if (cumulativeProb + prob <= 0.95)
            numInCi += 1.0;
        else if (cumulativeProb + prob > 0.95 && cumulativeProb <= 0.95)
            numInCi += (0.95 - cumulativeProb) / prob;

        cumulativeProb += prob;
        std::cout << std::setw(3) << i+1 << " (" << v[i].first << ") " << prob << " " << cumulativeProb << std::endl;
        }
        
    std::cout << "   * Credible Interval Size: " << numInCi << std::endl;
    calculateMeanPartition();
    pairProbabilities();
}

void ParameterSummaryPartition::print(Model* m) {

    std::string truePartition = getTruePartition(m);
        
    std::vector<std::pair<std::string, int> > v;
    for (auto& it : partitionsStringMap)
        v.push_back(it);
    sort(v.begin(), v.end(), cmp);

    int nameSize = (int)name.length() + 2;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "   * " << name << ": ";
    double cumulativeProb = 0.0;
    double numInCi = 0.0;
    bool foundTruePartitionInSample = false;
    for (int i=0; i<v.size(); i++)
        {
        if (i != 0)
            {
            std::cout << "   * ";
            for (int i=0; i<nameSize; i++)
                std::cout << " ";
            }
        double prob = (double)v[i].second / numSamples;
        if (cumulativeProb + prob <= 0.95)
            numInCi += 1.0;
        else if (cumulativeProb + prob > 0.95 && cumulativeProb <= 0.95)
            numInCi += (0.95 - cumulativeProb) / prob;

        std::cout << std::setw(3) << i+1 << " (" << v[i].first << ") " << prob << " " << cumulativeProb + prob;
        if (v[i].first == truePartition)
            {
            std::cout << " <- True Partition";
            bool inCi = false;
            if (cumulativeProb + prob <= 0.95)
                inCi = true;
            else if (cumulativeProb + prob > 0.95 && cumulativeProb <= 0.95)
                {
                double u = RandomVariable::randomVariableInstance().uniformRv();
                if (u < (0.95 - cumulativeProb) / prob)
                    inCi = true;
                }
            std::cout << ((inCi == true) ? " YES" : " NO");
            foundTruePartitionInSample = true;
            }
        std::cout << std::endl;
        
        cumulativeProb += prob;
        }
    if (foundTruePartitionInSample == false)
        {
        std::cout << "   * ";
        for (int i=0; i<nameSize; i++)
            std::cout << " ";
        std::cout << "Did not find true partition, " << truePartition << ", among sampled partitions" << std::endl;
        }
}

std::string ParameterSummaryPartition::summarize(void) {


    return "";
}
