// FFTMetaProgram.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <cmath>

constexpr auto M_PI = 3.14159265358979323846;   // pi;
using namespace std;

// Implement a taylor expansion of the sine function purely as a metaprogramming function
template<unsigned M, unsigned N, unsigned B, unsigned A>
struct SinCosSeries {
	static double value() {
		return 1 - (A*M_PI / B)*(A*M_PI / B) / M / (M + 1)*SinCosSeries<M + 2, N, B, A>::value();
	}
};

template<unsigned N, unsigned B, unsigned A>
struct SinCosSeries<N, N, B, A> {
	static double value() { return 1.; }
};

template<unsigned B, unsigned A, typename T = double>
struct Sin;

template<unsigned B, unsigned A>
struct Sin<B, A, float> {
	static float value() {
		return (A*M_PI / B)*SinCosSeries<2, 24, B, A>::value();
	}
};

template<unsigned B, unsigned A>
struct Sin<B, A, double> {
	static double value() {
		return (A*M_PI / B)*SinCosSeries<2, 34, B, A>::value();
	}
};

template<unsigned B, unsigned A, typename T = double>
struct Cos;

template<unsigned B, unsigned A>
struct Cos<B, A, float> {
	static float value() {
		return SinCosSeries<1, 23, B, A>::value();
	}
};

template<unsigned B, unsigned A>
struct Cos<B, A, double> {
	static double value() {
		return SinCosSeries<1, 33, B, A>::value();
	}
};

template<unsigned N, typename T = double>
class DanielsonLanczos {
	DanielsonLanczos<N / 2, T> next;
public:
	void apply(T* data) {
		next.apply(data);
		next.apply(data + N);
		// Computing only two sine functions and starting with 1, the next roots (wr,wi) are calculated recurrently from (wpr,wpi): (wr,wi)+=(wr,wi)*(wpr,wpi).
		T w_temp, temp_real, temp_imag, w_real, w_imag, w_prev_real, w_prev_imag;
		w_temp = -Sin<N, 1, T>::value(); // sin(M_PI / N);
		w_prev_real = -2 * w_temp*w_temp;
		w_prev_imag = -Sin<N, 2, T>::value(); // -sin(2 * M_PI / N) 
		w_real = 1.0;
		w_imag = 0.0;
		for (unsigned i = 0; i < N; i += 2) {


			temp_real = data[i + N] * w_real - data[i + N + 1] * w_imag;
			temp_imag = data[i + N] * w_imag + data[i + N + 1] * w_real;
			data[i + N] = data[i] - temp_real;
			data[i + N + 1] = data[i + 1] - temp_imag;
			data[i] += temp_real;
			data[i + 1] += temp_imag;

			// roots of unity recurrence relationship
			w_temp = w_real;
			w_real += w_real * w_prev_real - w_imag * w_prev_imag;
			w_imag += w_imag * w_prev_real + w_temp * w_prev_imag;
		}
	}
};

// trivial recursion case
template<typename T>
class DanielsonLanczos<1, T> {
public:
	void apply(T* data) {}
};

// N=4 terminal case
template<typename T>
class DanielsonLanczos<4, T> {
public:
	void apply(T* data) {
		T t_real = data[2];
		T t_imag = data[3];
		data[2] = data[0] - t_real;
		data[3] = data[1] - t_imag;
		data[0] += t_real;
		data[1] += t_imag;
		t_real = data[6];
		t_imag = data[7];
		data[6] = data[5] - t_imag;
		data[7] = t_real - data[4];
		data[4] += t_real;
		data[5] += t_imag;

		t_real = data[4];
		t_imag = data[5];
		data[4] = data[0] - t_real;
		data[5] = data[1] - t_imag;
		data[0] += t_real;
		data[1] += t_imag;
		t_real = data[6];
		t_imag = data[7];
		data[6] = data[2] - t_real;
		data[7] = data[3] - t_imag;
		data[2] += t_real;
		data[3] += t_imag;
	}
};

// N=2 specialized class
template<typename T>
class DanielsonLanczos<2, T> {
public:
	void apply(T* data) {
		T t_real = data[2];
		T t_imag = data[3];
		data[2] = data[0] - t_real;
		data[3] = data[1] - t_imag;
		data[0] += t_real;
		data[1] += t_imag;
	}
};
template<unsigned P, typename T = double>
class GFFT {
	enum { N = 1 << P };
	DanielsonLanczos<N, T> recursion;
	void scramble(T* data) {
		// reverse-binary reindexing
		int i, m, j = 1;
		for (i = 1; i < 2 * N; i += 2) {
			if (j > i) {
				swap(data[j - 1], data[i - 1]);
				swap(data[j], data[i]);
			}
			m = N;
			while (m >= 2 && j > m) {
				j -= m;
				m >>= 1;
			}
			j += m;
		}

	}
public:
	void fft(T* data) {
		scramble(data);
		recursion.apply(data);
	}
};

constexpr size_t log_2(size_t n)
{
	return ((n <= 2) ? 1 : 1 + log_2(n / 2));
}

int main()
{
	constexpr unsigned int N = 8; // number of elements in the signal	
	unsigned int array_len = 2 * N; // break the signal into real, complex pairs as adjacent array elements
	double* data = new double[array_len];
	// Real values of the signal
	data[0] = 1;
	data[2] = 1;
	data[4] = 2;
	data[6] = 2;
	data[8] = 1;
	data[10] = 1;
	data[12] = 0;
	data[14] = 0;

	// Only real numbers getting passed in
	for (unsigned int i = 0; i < N; i++) {
		data[2 * i + 1] = 0.0;
	}

	std::cout << "data: ";

	std::copy(data, data + array_len, std::ostream_iterator<double>(std::cout, " "));
	auto x = GFFT<log_2(N), double>();
	x.fft(data);
	std::cout << "\nFFT'd data:" << endl;
	for (unsigned int i = 0; i < N; ++i) {
		char sign = '+';
		if (data[2 * i + 1] < 0) sign = '-';
		cout << setw(10) << setprecision(5) << data[2 * i] << "\t" << sign
			<< abs(data[2 * i + 1]) << "i" << endl;
	}
}

