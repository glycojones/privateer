// Library for the YSBL program Privateer (PRogramatic Identification of Various Anomalies Toothsome Entities Experience in Refinement)
// Licence: LGPL - Please check Licence.txt for details.
//
// 2013-
// York Structural Biology Laboratory
// The University of York

// All credits for the implementation of these error functions into Privateer go to Marcin Wojdyr(https://github.com/wojdyr)

#ifndef PRIVATEER_ERROR_H_INCLUDED
#define PRIVATEER_ERROR_H_INCLUDED

#include <stdexcept>  // for runtime_error
#include <string>
#include <utility>    // for forward


[[noreturn]]
inline void fail(const std::string& msg) { throw std::runtime_error(msg); }

template<typename T, typename... Args> [[noreturn]]
void fail(std::string&& str, T&& arg1, Args&&... args) {
  str += arg1;
  fail(std::move(str), std::forward<Args>(args)...);
}
template<typename T, typename... Args> [[noreturn]]
void fail(const std::string& str, T&& arg1, Args&&... args) {
  fail(str + arg1, std::forward<Args>(args)...);
}


// unreachable() is used to silence GCC -Wreturn-type and hint the compiler
[[noreturn]] inline void unreachable() {
#if defined(__GNUC__) || defined(__clang__)
  __builtin_unreachable();
#elif defined(_MSC_VER)
  __assume(0);
#endif
}

#endif
