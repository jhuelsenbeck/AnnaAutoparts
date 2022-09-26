#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include "EigenSystem.hpp"
#include "Msg.hpp"



EigenSystem::EigenSystem(void) {

    activeVals = 0;

    ccIjk[0] = new std::complex<double>[64];
    ccIjk[1] = new std::complex<double>[64];
    pi[0].resize(4);
    pi[1].resize(4);
}

EigenSystem::EigenSystem(const EigenSystem& e) {

    activeVals = e.activeVals;
    for (int n=0; n<2; n++)
        {
        eigenValues[n] = e.eigenValues[n];
        pi[n] = e.pi[n];
        ccIjk[n] = new std::complex<double>[64];
        for (int i=0; i<64; i++)
            ccIjk[n][i] = e.ccIjk[n][i];
        }
}

EigenSystem::~EigenSystem(void) {

    delete [] ccIjk[0];
    delete [] ccIjk[1];
}

EigenSystem& EigenSystem::operator=(const EigenSystem& e) {

    if (this != &e)
        {
        activeVals = e.activeVals;
        for (int n=0; n<2; n++)
            {
            eigenValues[n] = e.eigenValues[n];
            pi[n] = e.pi[n];
            for (int i=0; i<64; i++)
                ccIjk[n][i] = e.ccIjk[n][i];
            }
        }
    return *this;
}

void EigenSystem::calculateEigenSystem(NucleotideSquareMatrix_t& Q) {

    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    Eigen::EigenSolver< NucleotideSquareMatrix_t > eigSolver;
    eigSolver.compute( Q, true );
    eigenValues[activeVals] = eigSolver.eigenvalues();
    Eigen::Matrix<std::complex<double>, 4, 4> eigenVectors = eigSolver.eigenvectors();
    Eigen::Matrix<std::complex<double>, 4, 4> inverseEigenVectors = eigenVectors.inverse();

    // calculate cc_ijk
    std::complex<double>* pc = ccIjk[activeVals];
    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                 *(pc++) = eigenVectors(i,k) * inverseEigenVectors(k,j);
}

std::vector<double> EigenSystem::calulateStationaryFrequencies(NucleotideSquareMatrix_t& Q) {
    
    Eigen::VectorXd f = Q.transpose().fullPivLu().kernel();
    f = f / f.sum();
    //std::cout << f << std::endl;
    
    std::vector<double> stationaryFrequencies(4);
    for (int i=0; i<4; i++)
        stationaryFrequencies[i] = f(i);
    
    double sum = sumVector(stationaryFrequencies);
    if ( fabs(1.0 - sum) > 0.001 )
        {
        std::cout << "sum = " << sum << std::endl;
        Msg::error("Stationary frequencies don't sum to one");
        }
    for (int i=0; i<4; i++)
        {
        if (f[i] < 0.0)
            {
            std::cout << f << std::endl;
            Msg::error("Negative stationary frequency");
            }
        }
    
    
    return stationaryFrequencies;
}

void EigenSystem::flipActiveValues(void) {

    if (activeVals == 0)
        activeVals = 1;
    else
        activeVals = 0;
}

std::complex<double>* EigenSystem::getCijk(void) {

    return ccIjk[activeVals];
}

Eigen::Matrix<std::complex<double>, 4, 1>& EigenSystem::getEigenValues(void) {

    return eigenValues[activeVals];
}

double EigenSystem::sumVector(std::vector<double>& v) {

    double sum = 0.0;
    for (int i=0; i<4; i++)
        sum += v[i];
    return sum;
}

void EigenSystem::testStationaryFrequencies(NucleotideSquareMatrix_t& Q) {

}

void EigenSystem::updateNonReversibleRateMatrix(std::vector<double>& rates) {

    NucleotideSquareMatrix_t Q;
    
    // fill in off diagonal components of rate matrix
    int k = 0;
    for (int i=0; i<4; i++)
        {
        for (int j=0; j<4; j++)
            {
            if (i != j)
                Q(i,j) = rates[k++];
            }
        }
    
    // fill in the diagonal elements of the rate matrix
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
    
    // rescale the rate matrix
    std::vector<double> f = calulateStationaryFrequencies(Q);
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        averageRate += -f[i] * Q(i,i);
    double scaleFactor = 1.0 / averageRate;
    Q *= scaleFactor;
    
    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    calculateEigenSystem(Q);
    setStationaryFrequencies(f);
}

void EigenSystem::updateReversibleRateMatrix(std::vector<double>& rates, std::vector<double>& f) {

    NucleotideSquareMatrix_t Q;

    // fill in the off diagonal components of the rate matrix
    for (int i=0, k=0; i<4; i++)
        {
        for (int j=i+1; j<4; j++)
            {
            Q(i,j) = rates[k  ] * f[j];
            Q(j,i) = rates[k++] * f[i];
            }
        }

    // fill in the diagonal elements of the rate matrix
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

    // rescale the rate matrix
    double averageRate = 0.0;
    for (int i=0; i<4; i++)
        averageRate += -f[i] * Q(i,i);
    double scaleFactor = 1.0 / averageRate;
    Q *= scaleFactor;
    
    // calculate the Eigenvalues and Eigenvectors and do some precomputation
    calculateEigenSystem(Q);
    setStationaryFrequencies(f);
}
