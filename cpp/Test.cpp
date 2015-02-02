#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <armadillo>
#include <pthread.h>
#include "relayHARQ.h"

int main()
{
    /*
    arma::umat B = relayHARQ::get_n_diff_bits(4);
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            std::cout << B(i, j) << std::endl;
        }
    }
	*/

    /*
    arma::vec x = arma::tan((arma::linspace(1, 10, 10) - 0.5) * M_PI / 10) ;
    for (int i = 0; i < 10; i++)
    {
        std::cout << x(i) << std::endl;
    }

    std::complex<double> j(0, 1);
    arma::cx_vec omega = 0.25 * (1 + x * std::complex<double>(0, 1));

    for (int i = 0; i < 10; i++)
    {
        std::cout << omega(i) << std::endl;
    }

    arma::cx_mat y(3, 3, arma::fill::zeros);
    y(0, 0) = pow(abs(conj(std::complex<double>(3, 4))), 2);
    std::cout << y(0, 0) << std::endl;
    */

    unsigned int Nbps = 1; // Number of bits per symbol, 1 for BPSK, 2 for QPSK, 3 for 8QAM, 4 for 16QAM and 5 for 32QAM
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


    /*
    unsigned int p = 26;
    unsigned int i = 27;
    unsigned int j = 28;
    unsigned int q = 29;
    unsigned int k = 30;
    unsigned int l = 31;
    unsigned long int idx = j * pow(Q, 5) + l * pow(Q, 4) + i * pow(Q, 3) + k * pow(Q, 2) + p * Q + q; // Should be enough to handle 32 QAM since it takes 30 bits to handle 32 ^ 6 - 1
    std::cout << idx << std::endl;
    std::vector<unsigned int> pijqkl = relayHARQ::idx2pijqkl(idx, Q, relayHARQ::order);
    for (int w = 0 ; w < 20; w++)
    {
        relayHARQ::next(pijqkl, Q, relayHARQ::order);
        std::cout << pijqkl[0] << pijqkl[1] << pijqkl[2] << pijqkl[3] << pijqkl[4] << pijqkl[5] <<  std::endl;
    }
    */

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

    /*
     * Compute PEP for the symbols
     */

    unsigned int N = 200;
    double xi = 0.25;

    unsigned long int idx_first = 0; // The index corresponding to the first pijqkl to compute
    unsigned long int idx_last = pow(Q, 6) - 1; // The index corresponding to the last pijqkl to compute
    std::cout << idx_last << std::endl;
    std::vector<unsigned int> pijqkl = relayHARQ::idx2pijqkl(idx_first, Q, relayHARQ::order); // Initialize the index [p, i, j, q, k, l]
    arma::umat B = relayHARQ::get_n_diff_bits(Nbps);
    std::vector<double> PEP_MGF_bit(idx_last - idx_first + 1, 0); // The vector that stores the PEP for bits in the right order


    for (unsigned long int idx = idx_first; idx < idx_last + 1; idx++)
    {
        if (pijqkl[0] != pijqkl[3] or pijqkl[1] != pijqkl[4] or pijqkl[2] != pijqkl[5])
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
        std::cout << pijqkl[0] <<  pijqkl[1] <<pijqkl[2] <<pijqkl[3] <<pijqkl[4] <<pijqkl[0] << std::endl;
        relayHARQ::next(pijqkl, Q, relayHARQ::order);

        std::cout << "idx = " << idx << ": "<< PEP_MGF_bit[idx - idx_first] << std::endl;
    }


    std::string filename = "q3aptest_" + std::to_string(Q) + "QAM_" + std::to_string((int)Eb2N0) + "dB.data";
    char * filename_cstr = new char [filename.length() + 1];
    std::strcpy(filename_cstr, filename.c_str());
    std::cout << filename_cstr << std::endl;
    std::ofstream datafile;
    datafile.open(filename_cstr, std::ios::out);
    datafile.precision(18);
    datafile.setf(std::ios::scientific);

    datafile.write(reinterpret_cast<char*>(&PEP_MGF_bit[0]), PEP_MGF_bit.size() * sizeof(double));
    datafile.close();
	return 0;
}
