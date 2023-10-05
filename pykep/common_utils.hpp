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

#include <kep3/detail/s11n.hpp>
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
        // c contains a user-defined plapythonic entity and either:
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

// Check if a mandatory method is present in a user-defined entity.
void check_mandatory_method(const py::object &o, const char *s, const char *target);

// Check if the user is trying to construct a pykep object from a type, rather than from an object.
// This is an easy error to commit, and it is sneaky because the callable_attribute() machinery in
// check_mandatory_method() will detect the methods of the *class* (rather than instance methods), and it will thus not
// error out.
void check_not_type(const py::object &o, const char *target);

// A simple wrapper for getters. It will try to:
// - get the attribute "name" from the object o,
// - call it without arguments,
// - extract an instance from the ret value and return it.
// If the attribute is not there or it is not callable, the value "def_value" will be returned.
template <typename RetType>
static RetType getter_wrapper(const py::object &o, const char *name, const RetType &def_value)
{
    auto a = callable_attribute(o, name);
    if (a.is_none()) {
        return def_value;
    }
    return py::cast<RetType>(a());
}

// Helpers to implement pickling on top of Boost.Serialization.
template <typename T>
inline py::tuple pickle_getstate_wrapper(const T &x)
{
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oa(oss);
        oa << x;
    }

    return py::make_tuple(py::bytes(oss.str()));
}

template <typename T>
inline T pickle_setstate_wrapper(const py::tuple &state)
{
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, ("The state tuple passed to the deserialization wrapper "
                                    "must have 1 element, but instead it has "
                                    + std::to_string(py::len(state)) + " element(s)")
                                       .c_str());
    }

    auto *ptr = PyBytes_AsString(state[0].ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "A bytes object is needed in the deserialization wrapper");
    }

    std::istringstream iss;
    iss.str(std::string(ptr, ptr + py::len(state[0])));
    T x;
    {
        boost::archive::binary_iarchive iarchive(iss);
        iarchive >> x;
    }

    return x;
}

} // namespace pykep

#endif