#include <iostream>
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

	return 0;
}
