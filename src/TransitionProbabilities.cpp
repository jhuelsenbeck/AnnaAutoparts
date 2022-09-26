#include <iomanip>
#include <iostream>
#include "Branch.hpp"
#include "Chunk.hpp"
#include "Model.hpp"
#include "Node.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"

#undef DEBUG_TIPROBS



TransitionProbabilities::TransitionProbabilities(Model* m, Chunk* c, int nb, int nc) {

    activeProbs = 0;
    model = m;
    myDataChunk = c;
    numBranches = nb;
    numGammaCategories = nc;
    eigens = new EigenSystem;
    
    for (int i=0; i<2; i++)
        {
        transitionProbabilities[i].resize(numGammaCategories);
        for (int k=0; k<numGammaCategories; k++)
            transitionProbabilities[i][k].resize(numBranches);
            
        for (int k=0; k<numGammaCategories; k++)
            for (int b=0; b<numBranches; b++)
                transitionProbabilities[i][k][b].setZero();
        }
}

TransitionProbabilities::~TransitionProbabilities(void) {

    delete eigens;
}

void TransitionProbabilities::flipActiveEigens(void) {

    eigens->flipActiveValues();
}

void TransitionProbabilities::flipActiveProbs(void) {

    if (activeProbs == 0)
        activeProbs = 1;
    else
        activeProbs = 0;
}

void TransitionProbabilities::print(void) {

    for (int i=0; i<transitionProbabilities[activeProbs].size(); i++)
        {
        for (int k=0; k<numGammaCategories; k++)
            {
            std::cout << "P(" << i << "," << k << "):" << std::endl;
            std::cout << transitionProbabilities[activeProbs][k][i] << std::endl;
            }
        }
}

void TransitionProbabilities::setTransitionProbabilities(Tree* t, std::vector<double> gammaCategoryRates) {

    Eigen::Matrix<std::complex<double>, 4, 1>& ceigenvalue = eigens->getEigenValues();
    std::complex<double>* ccIjk = eigens->getCijk();
    std::vector<std::complex<double> > ceigValExp(4);

    double treeLength = model->getTreeLength(myDataChunk->getId());
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        Node* pAnc = p->getAncestor();
        if (pAnc != NULL)
            {
            Branch* b = t->findBranch(p, pAnc);
            double v = b->getProportion() * treeLength;
            
            for (int k=0; k<numGammaCategories; k++)
                {
                double r = gammaCategoryRates[k];
                
                for (int s=0; s<4; s++)
                    ceigValExp[s] = exp(ceigenvalue[s] * v * r);

                std::complex<double>* ptr = ccIjk;
                for (int i=0; i<4; i++)
                    {
                    for (int j=0; j<4; j++)
                        {
                        std::complex<double> sum = std::complex<double>(0.0, 0.0);
                        for(int s=0; s<4; s++)
                            sum += (*ptr++) * ceigValExp[s];
                        transitionProbabilities[activeProbs][k][p->getIndex()](i,j) = (sum.real() < 0.0) ? 0.0 : sum.real();
                        }
                    }

                }
            
            }
        }

#   ifdef DEBUG_TIPROBS
    std::cout << std::fixed << std::setprecision(4);
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        Node* pAnc = p->getAncestor();
        
        if (pAnc != NULL)
            {
            Branch* b = t->findBranch(p, pAnc);
            double v = b->getLength();
            std::cout << "v = " << v << std::endl;

            for (int i=0; i<4; i++)
                {
                for (int k=0; k<numGammaCategories; k++)
                    {
                    for (int j=0; j<4; j++)
                        std::cout << transitionProbabilities[activeProbs][k][p->getIndex()](i,j) << " ";
                    std::cout << "| ";
                    }
                std::cout << std::endl;
                }
            }
            
        }
#   endif
}

void TransitionProbabilities::updateRateMatrix(std::vector<double>& rates, std::vector<double>& f) {

    eigens->updateReversibleRateMatrix(rates, f);
}

