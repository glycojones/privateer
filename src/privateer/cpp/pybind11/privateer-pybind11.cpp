// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York


#include <pybind11/pybind11.h>
#include <exception>      

namespace py=pybind11;
void init_ccp4mg(py::module &m);
void init_restraints(py::module &m);
void init_pyanalysis(py::module &m);


PYBIND11_MODULE(libprivateer, m) {
    m.doc() = "Python wrapper for libPrivateer(C++) exposed via pybind11.";

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const std::exception& e) {
            PyErr_SetString(PyExc_RuntimeError, e.what());
        }
    });

    init_ccp4mg(m);
    init_restraints(m);
    init_pyanalysis(m);
}

// Prepare a google doc describing the organisation of python cppmodule and how they would work in practice.
// How would the user import module, submodules.
// list cons and pros. 