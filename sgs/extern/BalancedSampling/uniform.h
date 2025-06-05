#ifndef STDUNIFORM_HEADER
#define STDUNIFORM_HEADER

//#include <uniform_int_distribution>
//#include <uniform_real_distribution>
#include <random>

//BAD CHANGE THIS IN FUTURE TO NOT USE GLOBAL VARS
static std::random_device rd;
static std::mt19937 gen(rd());

//asserting that on the current OS uniform_int_distribution can use size_t as a template param
//this should almost always be the case, but size_t is not required to be an unsigned int type,
//and using uniform_int_distribution where size-t is NOT an unsigned int type may yield
//undefined behavior.
static_assert(
	std::is_same_v<std::size_t, unsigned> || 
	std::is_same_v<std::size_t, unsigned long> || 
	std::is_same_v<std::size_t, unsigned long long>
);

inline double stduniform() {
 	std::uniform_real_distribution<double> dist(0.0, 1.0);
	return dist(gen);
};

inline double stduniform(double v) {
	std::uniform_real_distribution<double> dist(0.0, v);
	return dist(gen);
};

inline int intuniform(int N) {
	std::uniform_int_distribution<int> dist(0, N);
	return dist(gen);
};

inline size_t sizeuniform(size_t N) {
	std::uniform_int_distribution<size_t> dist(0, N);
  	return dist(gen);
};

#endif
