#ifndef Probability_hpp
#define Probability_hpp

#ifndef PI
#    define PI 3.141592653589793
#endif

#include <vector>
class RandomVariable;
class Tree;

static bool availableNormalRv = false;
static double extraNormalRv = 0.0;



namespace  Probability {

    namespace Beta {
    
        double  pdf(double alpha, double beta, double x);
        double  lnPdf(double alpha, double beta, double x);
        double  rv(RandomVariable* rng, double a, double b);
        double  cdf(double alpha, double beta, double x);
        double  quantile(double alpha, double beta, double x);
    }

    namespace ChiSquare {
    
        double  pdf(double v, double x);
        double  lnPdf(double v, double x);
        double  rv(RandomVariable* rng, double v);
        double  cdf(double v, double x);
        double  quantile(double prob, double v);
    }

    namespace Dirichlet {
    
        double  pdf(const std::vector<double> &a, const std::vector<double> &z);
        double  lnPdf(const std::vector<double> &a, const std::vector<double> &z);
        bool    rv(RandomVariable* rng, const std::vector<double> &a, std::vector<double> &z);
    }

    namespace Gamma {
    
        double  pdf(double alpha, double beta, double x);
        double  lnPdf(double alpha, double beta, double x);
        double  rv(RandomVariable* rng, double alpha, double beta);
        double  cdf(double alpha, double beta, double x);
        void    discretization(std::vector<double> &catRate, double a, double b, int nCats, bool median);
    }
    
    namespace Geometric {
    
        int     rv(RandomVariable* rng, double p);
    }

    namespace Exponential {
    
        double  pdf(double lambda, double x);
        double  lnPdf(double lambda, double x);
        double  rv(RandomVariable* rng, double lambda);
        double  cdf(double lambda, double x);
    }

    namespace Normal {
    
        double  pdf(double mu, double sigma, double x);
        double  lnPdf(double mu, double sigma, double x);
        double  rv(RandomVariable* rng) ;
        double  rv(RandomVariable* rng, double mu, double sigma) ;
        double  cdf(double mu, double sigma, double x);
        double  quantile(double mu, double sigma, double p);
    }
    
    namespace Uniform {
    
        double  pdf(double low, double high, double x);
        double  lnPdf(double low, double high, double x);
        double  rv(RandomVariable* rng, double low, double high);
        double  cdf(double low, double high, double x);
    }

    namespace Helper {
        
        double  beta(double a, double b);
        double  chebyshevEval(double x, const double *a, const int n);
        double  epsilon(void);
        double  gamma(double x);
        bool    isFinite(double x);
        double  lnFactorial(int n);
        double  lnBeta(double a, double b);
        double  lnGamma(double a);
        double  lnGammacor(double x);
        double  incompleteBeta(double a, double b, double x);
        double  incompleteGamma(double x, double alpha, double LnGamma_alpha);
        void    normalize(std::vector<double>& vec);
        double  pointNormal(double prob);
        double  rndGamma(RandomVariable* rng, double s, bool& err);
        double  rndGamma1(RandomVariable* rng, double s);
        double  rndGamma2(RandomVariable* rng, double s);
    }
}

#endif
