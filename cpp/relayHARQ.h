/**
* Functions and classes used to formulate the Q3AP problem in relay-HARQ communication system.
*
*
*  \author      Wenhao Wu
*  \version     0.1
*  \date        2015-01-30
 */

#ifndef RELAYHARQ_H_
#define RELAYHARQ_H_

#include <armadillo>
#include <complex>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

namespace relayHARQ {

    /**
     * Different QAM constellations
     */
    const std::vector<std::complex<double> > QAM2 {{-1, 0}, {1, 0}};
    const std::vector<std::complex<double> > QAM4 {{-0.707106781186548, 0.707106781186548}, {-0.707106781186548, -0.707106781186548}, {0.707106781186548, 0.707106781186548}, {0.707106781186548, -0.707106781186548}};
    const std::vector<std::complex<double> > QAM8 {{-1.22474487139159, 0.408248290463863}, {-1.22474487139159, -0.408248290463863}, {-0.408248290463863, 0.408248290463863}, {-0.408248290463863, -0.408248290463863}, {1.22474487139159, 0.408248290463863}, {1.22474487139159, -0.408248290463863}, {0.408248290463863, 0.408248290463863}, {0.408248290463863, -0.408248290463863}};
    const std::vector<std::complex<double> > QAM16 {{-0.948683298050514, 0.948683298050514}, {-0.948683298050514, 0.316227766016838}, {-0.948683298050514, -0.948683298050514}, {-0.948683298050514, -0.316227766016838}, {-0.316227766016838, 0.948683298050514}, {-0.316227766016838, 0.316227766016838}, {-0.316227766016838, -0.948683298050514}, {-0.316227766016838, -0.316227766016838}, {0.948683298050514, 0.948683298050514}, {0.948683298050514, 0.316227766016838}, {0.948683298050514, -0.948683298050514}, {0.948683298050514, -0.316227766016838}, {0.316227766016838, 0.948683298050514}, {0.316227766016838, 0.316227766016838}, {0.316227766016838, -0.948683298050514}, {0.316227766016838, -0.316227766016838}};
    const std::vector<std::complex<double> > QAM32 {{-0.670820393249937, 1.11803398874990}, {-0.223606797749979, 1.11803398874990}, {-0.670820393249937, -1.11803398874990}, {-0.223606797749979, -1.11803398874990}, {-1.11803398874990, 0.670820393249937}, {-1.11803398874990, 0.223606797749979}, {-1.11803398874990, -0.670820393249937}, {-1.11803398874990, -0.223606797749979}, {-0.223606797749979, 0.670820393249937}, {-0.223606797749979, 0.223606797749979}, {-0.223606797749979, -0.670820393249937}, {-0.223606797749979, -0.223606797749979}, {-0.670820393249937, 0.670820393249937}, {-0.670820393249937, 0.223606797749979}, {-0.670820393249937, -0.670820393249937}, {-0.670820393249937, -0.223606797749979}, {0.670820393249937, 1.11803398874990}, {0.223606797749979, 1.11803398874990}, {0.670820393249937, -1.11803398874990}, {0.223606797749979, -1.11803398874990}, {1.11803398874990, 0.670820393249937}, {1.11803398874990, 0.223606797749979}, {1.11803398874990, -0.670820393249937}, {1.11803398874990, -0.223606797749979}, {0.223606797749979, 0.670820393249937}, {0.223606797749979, 0.223606797749979}, {0.223606797749979, -0.670820393249937}, {0.223606797749979, -0.223606797749979}, {0.670820393249937, 0.670820393249937}, {0.670820393249937, 0.223606797749979}, {0.670820393249937, -0.670820393249937}, {0.670820393249937, -0.223606797749979}};

    /**
     * The order of pijqkl
     */
    const unsigned int order[6] = {3, 0, 4, 1, 5, 2};

    /**
     * Channel types
     */
    enum Channel {AWGN, RAYLEIGH, RICIAN};

    /**
     * The class that defines the a set of test parameter
     */
    class TestParam
    {
    public:
        Channel type;
        arma::cx_vec mu_h;
        arma::cx_vec sigma2_h;
        arma::cx_vec sigma2_eps;
        unsigned int Nbps;
        std::vector<std::complex<double> > constellation;
        double sigma2_v;
        unsigned int N;
        double xi;
        std::string saved_file;

        /**
         * The default constructor
         */
        TestParam(){};

        /*
         * The fundamental constructor
         */
        TestParam(Channel type, arma::cx_vec mu_h, arma::cx_vec sigma2_h, arma::cx_vec sigma2_eps, unsigned int Nbps, double Eb2N0, unsigned int N, double xi, std::string saved_file);

        /*
         * The simple perfect Rician channel constructor
         */
        TestParam(double K, std::complex<double> amp_RD2SD, unsigned int Nbps, double Eb2N0, unsigned int N, double xi, std::string saved_file);
    };

    /**
     * The function to
     */

    /**
     * Compute the moment generating function (MGF) for the quadaratic form of Gaussian random variables z'Az.
     */
    arma::cx_vec get_Psi_Gaussian(arma::cx_vec omega, arma::cx_vec mu, arma::cx_mat R, arma::cx_mat A);

    /**
     * Compute the pairwise error probability between 2 sets of symbols in the relay-HART model.
     */
    double get_PEP_symbol(arma::cx_vec x, arma::cx_vec y, arma::cx_vec mu_h, arma::cx_vec sigma2_h, arma::cx_vec sigma2_eps, double sigma2_v, int N, double xi);

    /**
     * Get the number of different bits between integers in [0, 2 ^ Nbps - 1]
     */
    arma::umat get_n_diff_bits(unsigned int Nbps);

    /**
     * Convert a index number (0 : Q ^ 6 - 1)  to the corresponding 6-number index [p, i, j, q, k, l]
     */
    std::vector<unsigned int> idx2pijqkl(unsigned long int idx, unsigned int Q, const unsigned int order[]);

    /**
     * Increase the 6-number index by 1.
     */
    void next(std::vector<unsigned int> & pijqkl, unsigned int Q, const unsigned int order[]);

}
#endif
