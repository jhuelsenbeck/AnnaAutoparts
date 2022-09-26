#include <iostream>
#include "RateMatrix.hpp"



RateMatrix::RateMatrix(int na) {

    numMatrices = na + 1;
    activeMatrix = 0;
    eigs.resize(numMatrices);
}

void RateMatrix::updateRateMatrix(std::vector<double>& r, std::vector<double>& f) {

    NucleotideSquareMatrix_t Q;

    // fill in the off-diagonal components of the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            Q(i,j) = r[k] * f[j];
            Q(j,i) = r[k] * f[i];
            k++;
            }
        }
        
    // fill in the diagonal components
    for (int i=0; i<4; i++)
        {
        double sum = 0.0;
        for (int j=0; j<4; j++)
            {
            if (i != j)
                sum += Q(i,j);
            }
        Q(i,i) = -sum;
        }
        
    // rescale the rate matrix so the average rate is 1.0
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        averageRate += -(f[i] * Q(i,i));
    double factor = 1.0 / averageRate;
    Q *= factor;
    
    std::cout << Q << std::endl;
    
    // calculate the eigen values and vectors of the rate matrix
    // and also set up the c_ijk vector for computing transition
    // probabilities
    eigs[activeMatrix].calculateEigenSystem(Q);
    eigs[activeMatrix].setStationaryFrequencies(f);
}
