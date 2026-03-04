#ifdef _WIN32
# ifdef LIBCCP4_EXPORT
#  define CCP4_IMEX __declspec(dllexport)
# else
#  define CCP4_IMEX __declspec(dllimport)
# endif
#else
# if (__GNUC__ > 4) || (__GNUC__ == 4 && (defined(__APPLE__) || __GNUC_MINOR__ >= 3))
#  define CCP4_IMEX __attribute__((visibility("default")))
# else
#  define CCP4_IMEX
# endif
#endif
