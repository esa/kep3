// Copyright 2023, 2024 Dario Izzo (dario.izzo@gmail.com), Francesco Biscani
// (bluescarni@gmail.com)
//
// This file is part of the kep3 library.
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef PYKEP_COMMON_UTILS_HPP
#define PYKEP_COMMON_UTILS_HPP

#include <sstream>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace pykep
{

namespace py = pybind11;

// Import and return the builtins module.
py::object builtins();

// Throw a Python exception.
[[noreturn]] void py_throw(PyObject *, const char *);

// Check if 'o' has a callable attribute (i.e., a method) named 's'. If so, it will
// return the attribute, otherwise it will return None.
py::object callable_attribute(const py::object &, const char *);

// Check if type is callable.
bool callable(const py::object &);

// Get the string representation of an object.
std::string str(const py::object &);

// Get the type of an object.
py::object type(const py::object &);

// Perform a deep copy of input object o.
py::object deepcopy(const py::object &);

// repr() via ostream.
template <typename T>
inline std::string ostream_repr(const T &x)
{
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

// Generic copy wrappers.
template <typename T>
inline T generic_copy_wrapper(const T &x)
{
    return x;
}

template <typename T>
inline T generic_deepcopy_wrapper(const T &x, const py::dict &)
{
    return x;
}

// Generic extract() wrappers.
template <typename C, typename T>
inline T *generic_cpp_extract(C &c, const T &)
{
    return c.template extract<T>();
}

template <typename C>
inline py::object generic_py_extract(C &c, const py::object &t)
{
    auto ptr = c.template extract<py::object>();
    if (ptr && (t.equal(type(*ptr)) || t.equal(builtins().attr("object")))) {
        // c contains a user-defined pythonic entity and either:
        // - the type passed in by the user is the exact type of the user-defined
        //   entity, or
        // - the user supplied as t the builtin 'object' type (which we use as a
        //   wildcard for any Python type).
        // Let's return the extracted object.
        return *ptr;
    }

    // Either the user-defined entity is not pythonic, or the user specified the
    // wrong type. Return None.
    return py::none();
}

} // namespace pykep

#endif