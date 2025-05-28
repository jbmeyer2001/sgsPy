#include <boost/python.hpp>
#include <Python.h>

char const* greet() {
	return "hello, world";
}

BOOST_PYTHON_MODULE(test_module) {
	using namespace boost::python;
	def("greet", greet);
}
