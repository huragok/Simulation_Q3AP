#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <armadillo>
#include <pthread.h>
#include <sys/time.h>
#include "relayHARQ.h"
#include <omp.h>

int main()
{
    struct timespec start, finish; // Time recording variables
    double elapsed;

    /**
     * Relay fading setup
     */
    arma::vec phi(4);
    phi(0) = 0, phi(1) = M_PI / 12, phi(2) = M_PI / 6, phi(3) = M_PI / 4; // Rician channel, the phase of the LOS component for the R-D link
    arma::vec amp(3);
    amp(0) = 1, amp(1) = sqrt(2), amp(2) = 2; // Rician channel, the amplitude of the LOS component for the R-D link
    arma::cx_vec relayLOS(12);
    for (int i_phi = 0; i_phi < 4; i_phi++)
    {
        for (int i_amp = 0; i_amp < 3; i_amp++)
        {
            relayLOS(i_phi * 3 + i_amp) = amp(i_amp) * std::exp(std::complex<double>(0, 1) * phi(i_phi));
        }
    }
    /**
     *  Simulation setup
     */
    double K = 10;
    unsigned int Nbps = 3; // Number of bits per symbol, 1 for BPSK, 2 for QPSK, 3 for 8QAM, 4 for 16QAM and 5 for 32QAM
    double Eb2N0 = 0; // Eb/N0 in dB
    unsigned int N0 = 200; // Size of the serial expansion
    double xi0 = 0.25;
    unsigned int Q0 = (unsigned int)pow(2, Nbps);

    std::vector<relayHARQ::TestParam> test_cases(12);
    for (int i_case = 0; i_case < 12; i_case++) // Generate all test cases.
    {
        std::string filename = "q3aptest_" + std::to_string(Q0) + "QAM_" + std::to_string((int)Eb2N0) + "dB_" + std::to_string((int)(K)) + "_case" + std::to_string(i_case) + ".data";
        test_cases[i_case] = relayHARQ::TestParam(K, relayLOS(i_case), Nbps, Eb2N0, N0, xi0, filename);
    }


    /*
     * Start generating the cost matrix based on the TestParam object
     */
    for (int i_case = 0; i_case < 12; i_case++)
    {
        arma::cx_vec mu_h = test_cases[i_case].mu_h;
        arma::cx_vec sigma2_h = test_cases[i_case].sigma2_h;
        arma::cx_vec sigma2_eps = test_cases[i_case].sigma2_eps;
        std::vector<std::complex<double> > constellation = test_cases[i_case].constellation;
        double sigma2_v = test_cases[i_case].sigma2_v;
        unsigned int N = test_cases[i_case].N;
        double xi = test_cases[i_case].xi;

        unsigned int Q = (unsigned int)pow(2, test_cases[i_case].Nbps);
        arma::umat B = relayHARQ::get_n_diff_bits(test_cases[i_case].Nbps);
        std::vector<double> PEP_MGF_bit(pow(Q, 6), 0); // The vector that stores the PEP for bits in the right order

        clock_gettime(CLOCK_MONOTONIC, &start);

        // Use OpenMP to parallelize, can use "num_threads()" to set the number of threads, or remove it to spawn the default number of threads
        #pragma omp parallel num_threads(8) firstprivate(mu_h, sigma2_h, sigma2_eps, constellation, sigma2_v, N, xi, Q, B) shared(PEP_MGF_bit)
        {
            #pragma omp single
            {
                int nthreads = omp_get_num_threads();
                std::cout << "Number of threads: " << nthreads << std::endl;
            }

            #pragma omp for
            for (unsigned long int idx = 0; idx < (int)pow(Q, 6); idx++)
            {
                std::vector<unsigned int> pijqkl = relayHARQ::idx2pijqkl(idx, Q, relayHARQ::order); // Initialize the index [p, i, j, q, k, l]
                if (pijqkl[0] != pijqkl[3] && pijqkl[1] != pijqkl[4] && pijqkl[2] != pijqkl[5])
                {
                    arma::cx_vec symbols_base(3, arma::fill::zeros);
                    symbols_base(0) = constellation[pijqkl[0]];
                    symbols_base(1) = constellation[pijqkl[1]];
                    symbols_base(2) = constellation[pijqkl[2]];

                    arma::cx_vec symbols_alt(3, arma::fill::zeros);
                    symbols_alt(0) = constellation[pijqkl[3]];
                    symbols_alt(1) = constellation[pijqkl[4]];
                    symbols_alt(2) = constellation[pijqkl[5]];

                    double PEP_MGF = relayHARQ::get_PEP_symbol(symbols_base, symbols_alt, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi);
                    PEP_MGF_bit[idx] = PEP_MGF * B(pijqkl[0], pijqkl[3]) / Q; // Although c++ vector is generally not thread safe, various sources have claimed that this operation should be okay. Not sure whether there is a racing problem.
                }

            }
        }
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        std::cout << "Case " << i_case << ": " << test_cases[i_case].saved_file << " generated." << std::endl;
        std::cout << "The LOS component of the R-D link is " << real(mu_h(2)) << " + " << imag(mu_h(2)) << "j" << std::endl;
        std::cout << "Elapsed time is: " << elapsed << "s" << std::endl;

        /**
         * Write to the file
         */
        char * filename_cstr = new char [test_cases[i_case].saved_file.length() + 1];
        std::strcpy(filename_cstr, test_cases[i_case].saved_file.c_str());
        std::ofstream datafile;
        datafile.open(filename_cstr, std::ios::out | std::ios::binary);

        datafile.setf(std::ios::scientific);
        datafile.precision(16);

        for (std::vector<double>::iterator itr(PEP_MGF_bit.begin()); itr != PEP_MGF_bit.end(); itr++)
        {
            datafile << *itr << "  ";
        }

        datafile.close();
    }

	return 0;
}
