/*
 *  \author    Wenhao Wu
 *  \version   0.1
 *  \date      2015-01-30  
 */

#include "relayHARQ.h"

namespace relayHARQ{

    arma::cx_vec get_Psi_Gaussian(arma::cx_vec omega, arma::cx_vec mu, arma::cx_mat R, arma::cx_mat A)
    {
        unsigned int N = omega.n_elem; // Size of serial expansion
        unsigned int d = mu.n_elem; // Size of the Gaussian random vector
        arma::cx_vec psi(N, arma::fill::zeros);
        arma::cx_mat eye_d(d, d, arma::fill::eye); // d-by-d identity matrix
        for (unsigned int n = 0; n < N; n++)
        {
            arma::cx_mat IRA = eye_d + omega(n) * R * A;
            arma::cx_vec IRA_mu = arma::solve(IRA, mu);
            arma::cx_mat temp = -omega(n) * mu.t() * A * IRA_mu;
            psi(n) = std::exp(temp(0)) / arma::det(IRA);
        }
        return psi;
    }

    double get_PEP_symbol(arma::cx_vec x, arma::cx_vec y, arma::cx_vec mu_h, arma::cx_vec sigma2_h, arma::cx_vec sigma2_eps, double sigma2_v, int N, double xi)
    {
        arma::vec tau = arma::tan((arma::linspace(1, N, N) - 0.5) * M_PI / N);
        arma::cx_vec omega = xi * (1 + tau * std::complex<double>(0, 1));

        // Phase 1
        arma::cx_vec mu0(3, arma::fill::zeros);
        mu0(0) = mu_h(0);

        arma::cx_mat R0(3, 3, arma::fill::zeros);
        R0(0, 0) = sigma2_h(0);
        R0(1, 1) = sigma2_eps(0);
        R0(2, 2) = sigma2_v;

        std::complex<double> e0 = x(0) - y(0);

        arma::cx_mat A0(3, 3, arma::fill::zeros);
        A0(0, 0) = pow(abs(e0), 2);
        A0(0, 1) = -conj(e0) * y(0);
        A0(0, 2) = conj(e0);
        A0(1, 0) = conj(A0(0, 1));
        A0(1, 1) = pow(abs(y(0)), 2) - pow(abs(x(0)), 2);
        A0(1, 2) = conj(e0);
        A0(2, 0) = conj(A0(0, 2));
        A0(2, 1) = conj(A0(1, 2));

        // Phase 2
        arma::cx_vec mu12(5, arma::fill::zeros);
        mu12(0) = mu_h(1);
        mu12(1) = mu_h(2);

        arma::cx_mat R12(5, 5, arma::fill::zeros);
        R12(0, 0) = sigma2_h(1);
        R12(1, 1) = sigma2_h(2);
        R12(2, 2) = sigma2_eps(1);
        R12(3, 3) = sigma2_eps(2);
        R12(4, 4) = sigma2_v;

        std::complex<double> e1 = x(1) - y(1);
        std::complex<double> e2 = x(2) - y(2);

        arma::cx_mat A12(5, 5, arma::fill::zeros);
        A12(0, 0) = pow(abs(e1), 2);
        A12(0, 1) = conj(e1) * e2;
        A12(0, 2) = -conj(e1) * y(1);
        A12(0, 3) = -conj(e1) * y(2);
        A12(0, 4) = conj(e1);
        A12(1, 0) = conj(A12(0, 1));
        A12(1, 1) = pow(abs(e2), 2);
        A12(1, 2) = -conj(e2) * y(1);
        A12(1, 3) = -conj(e2) * y(2);
        A12(1, 4) = conj(e2);
        A12(2, 0) = conj(A12(0, 2));
        A12(2, 1) = conj(A12(1, 2));
        A12(2, 2) = pow(abs(y(1)), 2) - pow(abs(x(1)), 2);
        A12(2, 3) = conj(y(1)) * y(2) - conj(x(1)) * x(2);
        A12(2, 4) = conj(e1);
        A12(3, 0) = conj(A12(0, 3));
        A12(3, 1) = conj(A12(1, 3));
        A12(3, 2) = conj(A12(2, 3));
        A12(3, 3) = pow(abs(y(2)), 2) - pow(abs(x(2)), 2);
        A12(3, 4) = conj(e2);
        A12(4, 0) = conj(A12(0, 4));
        A12(4, 1) = conj(A12(1, 4));
        A12(4, 2) = conj(A12(2, 4));
        A12(4, 3) = conj(A12(3, 4));

        arma::cx_vec psi = get_Psi_Gaussian(omega, mu0, R0, A0) % get_Psi_Gaussian(omega, mu12, R12, A12);

        double p = sum(real(psi) + tau % imag(psi)) / (2 * N);

        return p;
    }

    arma::umat get_n_diff_bits(unsigned int Nbps)
    {
        unsigned int Q = (unsigned int)(pow(2, Nbps));
        arma::umat B(Q, Q);
        for (unsigned int i = 0; i < Q; i++)
        {
            for (unsigned int j = 0; j < Q; j++)
            {
                unsigned int diff = i ^ j;

                for (B(i, j) = 0; diff; diff >>= 1)
                {
                    B(i, j) += diff & 1;
                }
            }
        }
        return B;
    }

    std::vector<unsigned int> idx2pijqkl(unsigned long int idx, unsigned int Q, const unsigned int order[])
    {
        std::vector<unsigned int> pijqkl(6, 0);

        for (int d = 0; d < 6; d++)
        {
            pijqkl[order[d]] = idx % Q;
            idx /= Q;
        }

        return pijqkl;
    }

    void next(std::vector<unsigned int> & pijqkl, unsigned int Q, const unsigned int order[])
    {
        unsigned int d = 0;
        while (pijqkl[order[d]] == Q - 1 and d < 6)
        {
            pijqkl[order[d]] = 0;
            d++;
        }
        if (d < 6)
        {
            pijqkl[order[d]]++;
        }

        return;
    }

}

