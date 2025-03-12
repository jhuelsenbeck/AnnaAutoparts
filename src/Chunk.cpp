#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Chunk.hpp"
#include "ConditionalLikelihoods.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include <unordered_map>

#undef NO_DATA


Chunk::Chunk(const Chunk& c) {

    alignment = c.alignment;
    model = c.model;
    id = c.id;
    numSites = c.numSites;
    numGammaCategories = c.numGammaCategories;

    condLikes = new ConditionalLikelihoods(alignment, numGammaCategories);
    transProbs = new TransitionProbabilities(model, this, 2*alignment->getNumTaxa()-1, numGammaCategories);
}

Chunk::Chunk(int myId, Alignment* a, Model* m, int nc) {

    alignment = a;
    model = m;
    id = myId;
    
    numSites = alignment->getNumSites();
    numGammaCategories = nc;
    
    condLikes = new ConditionalLikelihoods(alignment, numGammaCategories);
    transProbs = new TransitionProbabilities(model, this, 2*alignment->getNumTaxa()-1, numGammaCategories);
}

Chunk::~Chunk(void) {

    delete alignment;
    delete condLikes;
    delete transProbs;
}

int Chunk::getNumSitesForChunk(void) {

    return alignment->getNumSites();
}

void Chunk::flipActiveTransitionProbabilities(void) {

    transProbs->flipActiveProbs();
}

void Chunk::flipActiveEigens(void) {
    
    transProbs->flipActiveEigens();
}

double Chunk::lnLikelihood(void) {

#   if defined(NO_DATA)
    return 0.0;
#   endif

    // get the tree and a vector containing all of the transition probabilities
    Tree* t = model->getTree(id);
    std::vector<std::vector<NucleotideSquareMatrix_t> >& tp = transProbs->getTransitionProbabilities();
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    if (t->getRoot()->getDescendants().size() != 3)
        Msg::error("Expecting three descendants of the root");
    
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p->getIsLeaf() == false)
            {
            std::vector<Node*> desc = p->getDescendants();
            
            if (desc.size() == 2)
                {
                double* clL = condLikes->getConditionalLikelihoods(0, desc[0]->getIndex());
                double* clR = condLikes->getConditionalLikelihoods(0, desc[1]->getIndex());
                double* clP = condLikes->getConditionalLikelihoods(0, p->getIndex()      );
                for (int c=0; c<numSites; c++)
                    {
                    for (int k=0; k<numGammaCategories; k++)
                        {
                        NucleotideSquareMatrix_t& tpL = tp[k][desc[0]->getIndex()]; // check to see if this can be brought outside of the loop
                        NucleotideSquareMatrix_t& tpR = tp[k][desc[1]->getIndex()]; // over the number of sites...may be slow
                        for (int i=0; i<4; i++)
                            {
                            double sumL = 0.0, sumR = 0.0;
                            for (int j=0; j<4; j++)
                                {
                                sumL += tpL(i,j) * clL[j];
                                sumR += tpR(i,j) * clR[j];
                                }
                            clP[i] = sumL * sumR;
                            }
                        clP += 4;
                        clL += 4;
                        clR += 4;
                        }
                    }
                }
            else if (desc.size() == 3)
                {
                if (p != t->getRoot())
                    Msg::error("Only expecting a three-way split at the root of the tree");
                double* cl0 = condLikes->getConditionalLikelihoods(0, desc[0]->getIndex());
                double* cl1 = condLikes->getConditionalLikelihoods(0, desc[1]->getIndex());
                double* cl2 = condLikes->getConditionalLikelihoods(0, desc[2]->getIndex());
                double* clP = condLikes->getConditionalLikelihoods(0, p->getIndex()      );
                for (int c=0; c<numSites; c++)
                    {
                    for (int k=0; k<numGammaCategories; k++)
                        {
                        NucleotideSquareMatrix_t& tp0 = tp[k][desc[0]->getIndex()]; // check to see if this can be brought outside of the loop
                        NucleotideSquareMatrix_t& tp1 = tp[k][desc[1]->getIndex()]; // over the number of sites...may be slow
                        NucleotideSquareMatrix_t& tp2 = tp[k][desc[2]->getIndex()]; // over the number of sites...may be slow
                        for (int i=0; i<4; i++)
                            {
                            double sum0 = 0.0, sum1 = 0.0, sum2 = 0.0;
                            for (int j=0; j<4; j++)
                                {
                                sum0 += tp0(i,j) * cl0[j];
                                sum1 += tp1(i,j) * cl1[j];
                                sum2 += tp2(i,j) * cl2[j];
                                }
                            clP[i] = sum0 * sum1 * sum2;
                            }
                        clP += 4;
                        cl0 += 4;
                        cl1 += 4;
                        cl2 += 4;
                        }
                    }
                }
            else
                {
                Msg::error("Incorrect number of descendants!");
                }
            }
        }
        
    // calculate average probability at root
    std::vector<double> f = model->getBaseFrequencies(id);
    double lnL = 0.0;
    double catProb = 1.0 / numGammaCategories;
    double* clRoot = condLikes->getConditionalLikelihoods(0, t->getRoot()->getIndex());
    for (int c=0; c<numSites; c++)
        {
        double prob = 0.0;
        for (int k=0; k<numGammaCategories; k++)
            {
            for (int i=0; i<4; i++)
                prob += clRoot[i] * f[i] * catProb;
            clRoot += 4;
            }
        int numSitesOfPattern = alignment->getNumSitesOfPattern(c);
        lnL += numSitesOfPattern * log(prob);
        }
        
    return lnL;
}

void Chunk::print(void) {

    std::cout << "Data Subset ID: " << id << std::endl;
    alignment->print();
}

void Chunk::printTransitionProbabilities(void) {

    transProbs->print();
}

void Chunk::updateRateMatrix(void) {

    std::vector<double> freqs = model->getBaseFrequencies(id);
    std::vector<double> rates = model->getExchangeabilityRates(id);
    transProbs->updateRateMatrix(rates, freqs);
}

void Chunk::updateTransitionProbabilities(void) {

    Tree* t = model->getTree(id);
    std::vector<double> gammaRates = model->getGammaCategoryRates(id);
    transProbs->setTransitionProbabilities(t, gammaRates);
}
