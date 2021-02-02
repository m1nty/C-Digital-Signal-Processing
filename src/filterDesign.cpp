#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>

#define PI 3.14159265358979323846

// function for DFT (reused from previous experiment)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf)
{
	Xf.resize(x.size(), static_cast<std::complex<float>>(0, 0));
	for (auto m = 0; m < Xf.size(); m++) {
		for (auto k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to print a real vector (reused from previous experiment)
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (auto i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

// function to print a complex vector (reused from previous experiment)
void printComplexlVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (auto i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}


//***********************************TAKE HOME EXERCISE #2******************************

// Slice Function for sub-vectors
std::vector<float> slice(std::vector<float>arr,int low, int high)
{
    auto start = arr.begin() + low;
    auto end = arr.begin() + high + 1;
    std::vector<float>  result(high - low + 1);
    copy(start, end, result.begin());
    return result;
}

// function for DFT in segments
void dftSegments(const std::vector<float> &x, std::vector<std::complex<float>> &Xf_1024, std::vector<std::complex<float>> &Xf_256)
	{
			auto start_time_1024 = std::chrono::high_resolution_clock::now();
			DFT(x,Xf_1024);
			auto stop_time_1024 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> run_time_1024 = stop_time_1024-start_time_1024;
			std::cout << "1024 Point DFT ran for " << run_time_1024.count() << " milliseconds" << "\n";

			std::vector<std::complex<float>> bin_1;
			std::vector<std::complex<float>> bin_2;
			std::vector<std::complex<float>> bin_3;
			std::vector<std::complex<float>> bin_4;

			auto start_time_256= std::chrono::high_resolution_clock::now();
			DFT(slice(x, 0, 255), bin_1);
			DFT(slice(x, 256, 511), bin_2);
			DFT(slice(x, 512, 767), bin_3);
			DFT(slice(x, 768, 1023), bin_4);

			Xf_256.resize(256, static_cast<std::complex<float>>(0, 0));

			for (auto i = 0 ; i < 256 ; i++){
				Xf_256[i] = abs(((bin_1[i]) + (bin_2[i]) + (bin_3[i]) + (bin_4[i]) )/std::complex<float>(4, 0));
			}

			auto stop_time_256 = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double, std::milli> run_time_256 = stop_time_256-start_time_256;
			std::cout << "4 Segments of 256 Point DFT ran for " << run_time_256.count() << " milliseconds" << "\n";
	}

// function to generate a sine with N samples per second over interval
void generateSin(std::vector<float> &t, std::vector<float> &x, float Fs, float interval, float frequency = 7.0, float amplitude = 5.0, float phase = 0.0)
{
	// we do NOT allocate memory space explicitly
	// for the time (t) vector and sample (x) vector
	t.resize(0); x.resize(0);
	float dt = 1/Fs;
	for (auto i = 0.0; i < interval; i += dt) {
		// vector size increases when pushing new elements into it
		t.push_back(i);
		x.push_back(amplitude*std::sin(2*PI*frequency*i+phase));
	}
}

// function to mix an array of sines
void mixSin(const std::vector<std::vector<float>> &sv, std::vector<float> &mixed)
{
	// assumes at least one sine passed
	// assumes all input sines are of the same size
	for (auto i = 0.0; i < sv[0].size(); i ++) {
		float mixval = 0.0;
		// note: sv.size() returns the number of sines (or rows in 2D repr)
		// sv[0].size() returns the number of samples in a sine (or cols in 2D repr)
		for (auto k = 0; k < sv.size(); k++)
			mixval += sv[k][i];
		mixed.push_back(mixval);
	}
}


// function to record data in a format to be read by GNU plot
// the arguments are VERY specific to this usage in this experiment
// we have the time vector (t), a vector of sines (sv),
// input samples (x, i.e., mixed sines for this experiment),
// output samples (y, filtered input samples)
// frequency vectors (both Xf and Yf)
// note: the reference code does NOT do filtering, hence y and Yf are zeros by default
void plotMixedSinesSpectrum(const std::vector<float> &t, const std::vector<std::vector<float>> &sv, const std::vector<float> &x, const std::vector<float> &y, const std::vector<std::complex<float>> &Xf, const std::vector<std::complex<float>> &Yf)
{
	// write data in text format to be parsed by gnuplot
	const std::string filename = "../data/example.dat";
	std::fstream fd;  // file descriptor
	fd.open(filename, std::ios::out);
	fd << "#\tindex\tsine(0)\tsine(1)\tsine(2)\tdata in\tdata out\tspectrum in\tspectrum out\n";
	for (auto i = 0; i < t.size(); i++) {
		fd << "\t " << i << "\t";
		for (auto k = 0; k < sv.size(); k++)
			fd << std::fixed << std::setprecision(3) << sv[k][i] << "\t ";
		fd << x[i] << "\t "<< y[i] << "\t\t ";
		fd << std::abs(Xf[i])/Xf.size() << "\t\t " << std::abs(Yf[i])/Yf.size() <<"\n";
	}
	std::cout << "Generated " << filename << " to be used by gnuplot\n";
	fd.close();
}

// function to compute the impulse response "h" based on the sinc function
// see pseudocode from previous lab that was implemented in Python
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);
	float normCutoff = Fc/(Fs/2);

	for (auto i = 0 ; i<num_taps ; i++){
		if(i == (num_taps-1)/2){
			h[i] = normCutoff;
		} else{
			h[i] = normCutoff * (sin(PI*normCutoff*(i-(num_taps-1)/2)) / (PI*normCutoff*(i-(num_taps-1)/2)));
		}
		h[i] = h[i] * pow(sin(i*PI/num_taps),2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"; this is based on
// your Python code from the take-home exercise from the previous lab
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// allocate memory for the output (filtered) data
	y.resize(x.size()+h.size()-1, 0.0);
	int M = x.size();
	int N = h.size();

	for (auto n = 0 ; n < (M+N-1) ; n++){
		for (auto k = 0 ; k < (N+1) ; k++){

			if((n-k)>=0 && k<=(N-1) && (n-k)<=(M-1)){
				y[n] += x[n-k]*h[k];
			}
		}
		if (n==x.size()-1){
			break;
		}
	}
}

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


int main()
{

	float Fs = 1024.0;                    // samples per second
	float interval = 1.0;                 // number of seconds
	unsigned short int num_taps = 101;    // number of filter taps
	float Fc = 50.0;                      // cutoff frequency (in Hz)

	// declare a vector of vectors for multiple sines
	std::vector<std::vector<float>> sv;
	// declare time and sine vectors
	std::vector<float> t, sine;
	// note: there is no explicit memory allocation through vector resizing
	// vector memory space will increase via the push_back method

	// generate and store the first tone
	// check the function to understand the order of arguments
	generateSin(t, sine, Fs, interval, 40.0, 5.0, 0.0);
	sv.push_back(sine);
	// generate and store the second tone
	generateSin(t, sine, Fs, interval, 100.0, 2.0, 0.0);
	sv.push_back(sine);
	// generate and store the third tone
	generateSin(t, sine, Fs, interval, 120.0, 3.0, 0.0);
	sv.push_back(sine);

	// declare the mixed sine vector and mix the three tones
	std::vector<float> x;
	mixSin(sv, x);
	// // printRealVector(x);
	//
	// // declare a vector of complex values for DFT; no memory is allocated for it
	// std::vector<std::complex<float>> Xf;
	// DFT(x, Xf);
	// // printComplexlVector(Xf);
	//
	// // generate the impulse response h
	// // convolve it with the input data x
	// // in order to produce the output data y
	// std::vector<float> h;              // impulse response
	// impulseResponseLPF(Fs, Fc, num_taps, h);
	// std::vector<float> y;              // filter out
	// convolveFIR(y, x, h);
	//
	// // compute DFT of the filtered data
	// std::vector<std::complex<float>> Yf;
	// DFT(y, Yf);
	//
	// // prepare the data for gnuplot
	//plotMixedSinesSpectrum(t, sv, x, y, Xf, Yf);
	//
	// // naturally, you can comment the line below once you are comfortable to run gnuplot
	// std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";



	//***********************************TAKE HOME EXERCISE #2******************************

	std::vector<std::complex<float>> DFT_res_1024;
	std::vector<std::complex<float>> DFT_res_256;

	std::vector<float> h2;              // impulse response
	impulseResponseLPF(Fs, Fc, num_taps, h2);
	std::vector<float> y2;              // filter out
	convolveFIR(y2, x, h2);

	dftSegments(y2, DFT_res_1024, DFT_res_256);

	// prepare the data for gnuplot
	plotMixedSinesSpectrum(t, sv, x, y2, DFT_res_256, DFT_res_1024);


	return 0;
}
