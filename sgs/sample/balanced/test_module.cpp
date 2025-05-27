#include <iostream>
#include <Python.h>

static PyObject* foo(PyObject* self, PyObject* dataset) {
	std::cout << dataset << std::endl;
	return 0;
}

static PyMethodDef methods[] = {
	{"test", (PyCFunction)foo, METH_O, NULL}
};

static struct PyModuleDef module = {
	PyModuleDef_HEAD_INIT,
	"test_module",
	NULL,
	-1,
	methods,
};

PyMODINIT_FUNC PyInit_test_module(void) {
	return PyModule_Create(&module);
}
