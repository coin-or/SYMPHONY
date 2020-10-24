/*
 * Include file for the configuration of SYMPHONY.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN packages or third party code) are set
 * here.  The project maintainer needs to remember to update this file
 * and choose reasonable defines.  A user can modify the default
 * setting by editing this file here.
 *
 */

#ifndef __SYMPHONY_CONFIG_H__
#define __SYMPHONY_CONFIG_H__

#ifdef HAVE_CONFIG_H
#ifdef SYMPHONYLIB_BUILD
#include "config.h"

/* overwrite SYMPHONYLIB_EXPORT from config.h when building SYMPHONY
 * we want it to be __declspec(dllexport) when building a DLL on Windows
 * we want it to be __attribute__((__visibility__("default"))) when building
 * with GCC, so user can compile with -fvisibility=hidden
 */
#ifdef DLL_EXPORT
#undef SYMPHONYLIB_EXPORT
#define SYMPHONYLIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
#undef SYMPHONYLIB_EXPORT
#define SYMPHONYLIB_EXPORT __attribute__((__visibility__("default")))
#endif

#else
#include "config_sym.h"
#endif

#else /* HAVE_CONFIG_H */

#ifdef SYMPHONYLIB_BUILD
#include "config_default.h"
#else
#include "config_sym_default.h"
#endif

#endif /* HAVE_CONFIG_H */

#endif /*__SYMPHONY_CONFIG_H__*/
