#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>

#define PI 3.14159265358979323846

// // input x is qualified as const in order not to be changed within the function
// although C++ is case sensitive it is "safer" to use Xf instead of X in order to avoid confusion
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	// allocate space for the vector holding the frequency bins
	// initialize all the values in the frequency vector to complex 0
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	// "auto" keyword is used for type deduction (or inference) in C++
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

void IDFT(const std::vector<std::complex<float>> &X, std::vector<std::complex<float>> &xt) {
	// allocate space for the vector holding the frequency bins
	// initialize all the values in the frequency vector to complex 0
	// xt.resize(X.size(), static_cast<std::vector<float>>(0));
	xt.resize(X.size(), static_cast<std::complex<float>>(0, 0));
	// "auto" meyword is used for type deduction (or inference) in C++
	for (auto k = 0; k < xt.size(); k++) {
		for (auto m = 0; m < X.size(); m++) {
				std::complex<float> expval(0, 2*PI*(m*k) / X.size());
				xt[k] += X[m] * std::exp(expval);
		}
		xt[k] /= X.size();
	}
}

// function to generate N random float values
// x is the input/output vector
// N is the number of samples
// random values are between -max and max
// precision should be capped to 3 (in the context of this type of experiments)
void generateRandomSamples(std::vector<float> &x, unsigned int N, unsigned short int max, unsigned char precision)
{
	// allocate space for the vector holding the frequency bins
	x.resize(N);
	int int_radom_max = 2*(max * static_cast<int>(pow(10,precision)));
	for (auto i = 0; i < x.size(); i++) {
		// static casting is generally safer and better for both
		// type checking as well as for documentation purposes
		x[i] = (static_cast<float>(std::rand() % int_radom_max));
		x[i] = (x[i] - (int_radom_max/2))/pow(10,precision);
	}
}

// function to print a real vector
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (auto i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";

	// starting from C++11 you can write
	// for (auto elem : x)
	//   std::cout << elem << " ";

	std::cout << "\n";
}

// function to print a complex vector
void printComplexlVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (auto i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

int main()
{
	std::cout << "In-Lab Experiments\n" << '\n';

	// use initialize the seed for the random number generator using current time
	int seed = std::time(0x0);
	std::cout << "Starting from seed " << seed << "\n";

	// declare a vector of real values; no memory is allocated at this time
	std::vector<float> x;
	// generate 32 samples between -10 and 10; for extra flexibility
	// the last argument gives precision in terms of fraction digits
	// note: memory for x is allocated within the function called below
	generateRandomSamples(x, 128, 10, 2);
	// print a vector of real numbers
	// printRealVector(x);

	// declare a vector of complex values; no memory is allocated for it
	std::vector<std::complex<float>> Xf;
	std::vector<std::complex<float>> xt;

	// perform DFT of x to produce Xf
	// we measure the execution time using the "chrono" class
	auto start_time = std::chrono::high_resolution_clock::now();
 	DFT(x, Xf);

	auto stop_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> DFT_run_time = stop_time-start_time;
	// print a vector of complex numbers
	// printComplexlVector(Xf);
	std::cout << "DFT ran for " << DFT_run_time.count() << " milliseconds" << "\n";

	//In-Lab Experiment #1
	auto start_time2 = std::chrono::high_resolution_clock::now();
	IDFT(Xf, xt);
	auto stop_time2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> IDFT_run_time = stop_time2-start_time2;
	std::cout << "IDFT ran for " << IDFT_run_time.count() << " milliseconds" << "\n";



	//***********************************TAKE HOME EXERCISE #1******************************

	std::cout << "\nTakehome Exercise #1\n" << '\n';

  std::vector<float> samples;
	std::vector<std::complex<float>> Xf_res;
	std::vector<std::complex<float>> xt_res;
	for (auto N = 128; N <= pow(2,14) ; N*=2){
		generateRandomSamples(samples, N, 10, 2);
		auto start_times = std::chrono::high_resolution_clock::now();
		DFT(samples, Xf_res);
		IDFT(Xf,xt_res);
		auto stop_times = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double, std::milli> Run_Times = stop_times-start_times;
		std::cout << "DFT/IDFT for N = " << N << " ran for: " << Run_Times.count() << " milliseconds" << "\n";
	}

	return 0;
}
