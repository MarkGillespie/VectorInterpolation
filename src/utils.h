#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/utilities/vector3.h"

#include <iostream>
#include <sstream>
#include <string>

using geometrycentral::Vector3;
using std::string;

// cotan
double cot(double x);

// Custom runtime error which will report the offending file and line number
class verbose_runtime_error : public std::runtime_error {
    std::string msg;

  public:
    verbose_runtime_error(const std::string& arg, const char* file, int line);
    ~verbose_runtime_error() throw() {}
    const char* what() const throw() { return msg.c_str(); }
};

#define throw_verbose_runtime_error(arg)                                       \
    throw verbose_runtime_error(arg, __FILE__, __LINE__);

#define verbose_assert(arg, msg)                                               \
    if (!(arg)) throw_verbose_runtime_error(msg);

//== Linear algebra helpers
Eigen::Vector3d toEigen(const Vector3& v);
Vector3 fromEigen(const Eigen::Vector3d& v);
Eigen::Matrix3d outer(const Vector3& v, const Vector3& w); // Returns vw^T

//== Convenient printing helpers

string getFilename(string filePath, bool withExtension = true,
                   char seperator = '/');

// Verbose endl
#define vendl                                                                  \
    "\t\t(" << getFilename(__FILE__) << ":" << __LINE__ << ")" << std::endl

// == Horrible macro to print variable name + variable
// https://stackoverflow.com/a/6623090
#define WATCHSTR_WNAME(os, name, a)                                            \
    do {                                                                       \
        (os) << (name) << " is value " << (a) << vendl;                        \
    } while (false)

#define WATCHSTR(os, a) WATCHSTR_WNAME((os), #a, (a))
#define WATCH(a) WATCHSTR_WNAME(std::cout, #a, (a))
#define WATCH2(a, b)                                                           \
    WATCH(a);                                                                  \
    WATCH(b)
#define WATCH3(a, b, c)                                                        \
    WATCH(a);                                                                  \
    WATCH2(b, c)
#define WATCH4(a, b, c, d)                                                     \
    WATCH(a);                                                                  \
    WATCH3(b, c, d)
#define WATCH5(a, b, c, d, e)                                                  \
    WATCH(a);                                                                  \
    WATCH4(b, c, d, e)

#define HERE()                                                                 \
    do {                                                                       \
        std::cout << "HERE at " << __LINE__ << " in " << getFilename(__FILE__) \
                  << ":" << __FUNCTION__ << std::endl;                         \
    } while (false)
#define HERE1(name)                                                            \
    do {                                                                       \
        std::cout << "HERE (" << (name) << ") at " << __LINE__ << " in "       \
                  << getFilename(__FILE__) << ":" << __FUNCTION__              \
                  << std::endl;                                                \
    } while (false)

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data);

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& data);

template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& data);

template <typename T>
std::ostream& operator<<(std::ostream& out, const Eigen::Triplet<T>& data);


#include "utils.ipp"
