#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <armadillo>
#include <pthread.h>
#include <sys/time.h>
#include "relayHARQ.h"

int main()
{
    struct timespec start, finish;
    double elapsed;

    unsigned int Nbps = 2; // Number of bits per symbol, 1 for BPSK, 2 for QPSK, 3 for 8QAM, 4 for 16QAM and 5 for 32QAM
    unsigned int Q = (unsigned int)pow(2, Nbps);
    std::vector<std::complex<double> > constellation;
    switch (Nbps)
    {
        case 1:
            constellation = relayHARQ::QAM2;
            break;
        case 2:
            constellation = relayHARQ::QAM4;
            break;
        case 3:
            constellation = relayHARQ::QAM8;
            break;
        case 4:
            constellation = relayHARQ::QAM16;
            break;
        case 5:
            constellation = relayHARQ::QAM32;
            break;
        default:
            std::cerr << "Nbps set incorrectly!" << std::endl;
            exit(1);
	}

    /**
     * Channel settings (Let us try AWGN and Rayleigh fading first)
     */

    arma::cx_vec mu_h(3);
    arma::cx_vec sigma2_h(3);
    arma::cx_vec sigma2_eps(3);
    relayHARQ::Channel channel = relayHARQ::AWGN;
    switch (channel)
    {
        case relayHARQ::AWGN:
            mu_h(0) = 1, mu_h(1) = 1, mu_h(2) = 1;
            sigma2_h(0) = 0, sigma2_h(1) = 0, sigma2_h(2) = 0;
            sigma2_eps(0) = 0, sigma2_eps(1) = 0, sigma2_eps(2) = 0;
            break;
        case relayHARQ::RAYLEIGH:
            mu_h(0) = 0, mu_h(1) = 0, mu_h(2) = 0;
            sigma2_h(0) = 1, sigma2_h(1) = 1, sigma2_h(2) = 1;
            sigma2_eps(0) = 0, sigma2_eps(1) = 0, sigma2_eps(2) = 0;
            break;
        default:
            std::cerr << "Unknown channel type!" << std::endl;
            exit(1);
    }

    double Eb2N0 = 0; // Eb/N0 in dB
    double sigma2_v = 1.0 / (Nbps * pow(10, Eb2N0 / 10)); // The noise variance


    // Compute PEP for the symbols


    unsigned int N = 200;
    double xi = 0.25;

    unsigned long int idx_first = 0; // The index corresponding to the first pijqkl to compute
    unsigned long int idx_last = pow(Q, 6) - 1; // The index corresponding to the last pijqkl to compute
    std::vector<unsigned int> pijqkl = relayHARQ::idx2pijqkl(idx_first, Q, relayHARQ::order); // Initialize the index [p, i, j, q, k, l]
    arma::umat B = relayHARQ::get_n_diff_bits(Nbps);
    std::vector<double> PEP_MGF_bit(idx_last - idx_first + 1, 0); // The vector that stores the PEP for bits in the right order

    clock_gettime(CLOCK_MONOTONIC, &start);
    for (unsigned long int idx = idx_first; idx < idx_last + 1; idx++)
    {
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
            PEP_MGF_bit[idx - idx_first] = PEP_MGF * B(pijqkl[0], pijqkl[3]) / Q;
        }
        relayHARQ::next(pijqkl, Q, relayHARQ::order);

        //std::cout << "idx = " << idx << ": "<< PEP_MGF_bit[idx - idx_first] << std::endl;
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    std::cout << "Elapsed time is: " << elapsed << "s" << std::endl;

    std::string filename = "q3aptest_" + std::to_string(Q) + "QAM_" + std::to_string((int)Eb2N0) + "dB.data";
    char * filename_cstr = new char [filename.length() + 1];
    std::strcpy(filename_cstr, filename.c_str());
    //std::cout << filename_cstr << std::endl;
    std::ofstream datafile;
    datafile.open(filename_cstr, std::ios::out | std::ios::binary);

    datafile.setf(std::ios::scientific);
    datafile.precision(16);

    for (std::vector<double>::iterator itr(PEP_MGF_bit.begin()); itr != PEP_MGF_bit.end(); itr++)
    {
        datafile << *itr << "  ";
    }

    datafile.close();

	return 0;
}
