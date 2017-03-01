# -----------------------------------------------------------------
# $Revision: 1.9 $
# $Date: 2009/03/25 23:10:50 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for FCVODE module
#
# src/cvode/fcmix/Makefile.  Generated from Makefile.in by configure.
# -----------------------------------------------------------------

SHELL = /bin/bash

srcdir       = .
builddir     = .
abs_builddir = /home/brandon/GenFit/NFsim_v1.11/Network3/cvode-2.6.0/src/cvode/fcmix
top_builddir = ../../..
prefix       = /home/brandon/GenFit/NFsim_v1.11/Network3
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

INSTALL        = /usr/bin/install -c
INSTALL_LIB    = ${INSTALL}
INSTALL_HEADER = ${INSTALL} -m 644

LIBTOOL      = $(SHELL) $(top_builddir)/libtool
LIBTOOL_DEPS = config/ltmain.sh

CPP      = cc -E
CPPFLAGS = 
CC       = cc
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = -lm 

LAPACK_ENABLED = no

top_srcdir = $(srcdir)/../../..

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/src/cvode -I$(top_builddir)/include

LIB_REVISION = 0:1:0

FCVODE_LIB    = libsundials_fcvode.la

FCVODE_SRC_FILES = fcvode.c fcvband.c fcvdense.c fcvjtimes.c fcvpreco.c fcvbbd.c fcvbp.c fcvroot.c fcvewt.c
FCVODE_BL_SRC_FILES = fcvlapack.c fcvlapband.c fcvlapdense.c 

FCVODE_OBJ_FILES = $(FCVODE_SRC_FILES:.c=.o)
FCVODE_BL_OBJ_FILES = $(FCVODE_BL_SRC_FILES:.c=.o)

FCVODE_LIB_FILES = $(FCVODE_SRC_FILES:.c=.lo)
FCVODE_BL_LIB_FILES = $(FCVODE_BL_SRC_FILES:.c=.lo)

mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs

# ----------------------------------------------------------------------------------------------------------------------

all: $(FCVODE_LIB)

$(FCVODE_LIB): $(FCVODE_LIB_FILES)
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	  make lib_with_bl;                          \
	else                                         \
	  make lib_without_bl;                       \
	fi

lib_without_bl: $(FCVODE_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(FCVODE_LIB) $(FCVODE_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -static -version-info $(LIB_REVISION)

lib_with_bl: $(FCVODE_LIB_FILES) $(FCVODE_BL_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(FCVODE_LIB) $(FCVODE_LIB_FILES) $(FCVODE_BL_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -static -version-info $(LIB_REVISION)

install: $(FCVODE_LIB)
	$(mkinstalldirs) $(libdir)
	$(LIBTOOL) --mode=install $(INSTALL_LIB) $(FCVODE_LIB) $(libdir)

uninstall:
	$(LIBTOOL) --mode=uninstall rm -f $(libdir)/$(FCVODE_LIB)

clean:
	$(LIBTOOL) --mode=clean rm -f $(FCVODE_LIB)
	rm -f $(FCVODE_LIB_FILES)
	rm -f $(FCVODE_BL_LIB_FILES)
	rm -f $(FCVODE_OBJ_FILES)
	rm -f $(FCVODE_BL_OBJ_FILES)

distclean: clean
	rm -f Makefile

fcvode.lo: $(srcdir)/fcvode.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvode.c
fcvewt.lo: $(srcdir)/fcvewt.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvewt.c
fcvband.lo: $(srcdir)/fcvband.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvband.c
fcvdense.lo: $(srcdir)/fcvdense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvdense.c
fcvlapack.lo: $(srcdir)/fcvlapack.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvlapack.c
fcvlapband.lo: $(srcdir)/fcvlapband.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvlapband.c
fcvlapdense.lo: $(srcdir)/fcvlapdense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvlapdense.c
fcvjtimes.lo: $(srcdir)/fcvjtimes.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvjtimes.c
fcvpreco.lo: $(srcdir)/fcvpreco.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvpreco.c
fcvbbd.lo: $(srcdir)/fcvbbd.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvbbd.c
fcvbp.lo: $(srcdir)/fcvbp.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvbp.c
fcvroot.lo: $(srcdir)/fcvroot.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fcvroot.c

libtool: $(top_builddir)/$(LIBTOOL_DEPS)
	@cd ${top_builddir} ;                 \
	 ${SHELL} ./config.status --recheck ; \
	 cd ${abs_builddir}
