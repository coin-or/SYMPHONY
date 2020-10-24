/* Copyright (C) 2020
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * $Id: OsiClpConfig.h 2568 2020-01-08 08:57:16Z stefan $
 *
 * Include file for the configuration of OsiClp.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header file is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN-OR packages or third party code) are set
 * by the files config_*default.h. The project maintainer needs to remember
 * to update these files and choose reasonable defines.
 * A user can modify the default setting by editing the config_*default.h
 * files.
 */

#ifndef __OSISYMCONFIG_H__
#define __OSISYMCONFIG_H__

#ifdef HAVE_CONFIG_H
#ifdef OSISYMPHONYLIB_BUILD

#ifdef DLL_EXPORT
# define OSISYMPHONYLIB_EXPORT __declspec(dllexport)
#elif defined(__GNUC__) && __GNUC__ >= 4
# define OSISYMPHONYLIB_EXPORT __attribute__((__visibility__("default")))
#else
# define OSISYMPHONYLIB_EXPORT
#endif

#else
#include "config_osisym.h"
#endif

#else /* HAVE_CONFIG_H */

#ifndef OSISYMPHONYLIB_EXPORT
# if defined(_WIN32) && defined(DLL_EXPORT)
#  ifdef OSISYMPHONYLIB_BUILD
#   define OSISYMPHONYLIB_EXPORT __declspec(dllexport)
#  else
#   define OSISYMPHONYLIB_EXPORT __declspec(dllimport)
#  endif
# elif defined(__GNUC__) && __GNUC__ >= 4
#  define OSISYMPHONYLIB_EXPORT __attribute__((__visibility__("default")))
# else
#  define OSISYMPHONYLIB_EXPORT
# endif
#endif

#endif /* HAVE_CONFIG_H */

#endif /*__OSISYMCONFIG_H__*/
