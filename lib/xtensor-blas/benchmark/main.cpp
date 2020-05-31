// /***************************************************************************
// * Copyright (c) 2016, Johan Mabille, Sylvain Corlay and Wolf Vollprecht    *
// *                                                                          *
// * Distributed under the terms of the BSD 3-Clause License.                 *
// *                                                                          *
// * The full license is in the file LICENSE, distributed with this software. *
// ****************************************************************************/

#include <iostream>
#include <benchmark/benchmark.h>

#include "benchmark_blas.hpp"

#ifdef WITH_OPENBLAS
void blas_stats()
{
	std::cout << "\nSTATS FOR OPENBLAS\n---------------------------------------------------------------\n\n";

	openblas_set_num_threads(4);
	std::cout << "NUMBER OF THREADS:          " << openblas_get_num_threads() << std::endl;

	// Get the number of physical processors (cores)
	std::cout << "NUMBER OF PROCESSORS:       " << openblas_get_num_procs() << std::endl;
	// Get the build configure on runtime.
	std::cout << "CONFIG:                     " << openblas_get_config() << std::endl;

	/*Get the CPU corename on runtime.*/
	std::cout << "CORE NAME:                  " << openblas_get_corename() << std::endl;

	// Get the parallelization type which is used by OpenBLAS
	std::cout << "PARALLEL:                   " << openblas_get_parallel() << std::endl;
	std::cout << "\n\n";
}
#else
void blas_stats(){};
#endif

// Custom main function to print BLAS config
int main(int argc, char** argv)
{
	blas_stats();
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) return 1;
    benchmark::RunSpecifiedBenchmarks();
}