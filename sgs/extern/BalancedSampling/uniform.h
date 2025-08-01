#ifndef STDUNIFORM_HEADER
#define STDUNIFORM_HEADER

#include <random>

//BAD CHANGE THIS IN FUTURE TO NOT USE GLOBAL VARS
static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution dist(0.0, 1.0);

inline double stduniform() {
	return dist(gen);
}

inline double stduniform(double v) {
	return dist(gen) * v;
};

inline int intuniform(int N) {
	return (N == 0 || N == 1) ? 0 : (int)stduniform((double)N);
};

inline size_t sizeuniform(size_t N) {
	return (N == 0 || N == 1) ? 0 : (size_t)stduniform((double)N);
};

#endif
