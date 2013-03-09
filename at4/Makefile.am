## Process this file with automake to produce Makefile.in
# Note: Makefile is automatically generated from Makefile.in by the configure
# script, and Makefile.in is generated from Makefile.am by automake.

AM_CPPFLAGS = -I$(top_srcdir)/include -DGMXLIBDIR=\"$(datadir)/top\"


noinst_LTLIBRARIES  = libgmxpreprocess@LIBSUFFIX@.la

LDADD = libgmxpreprocess@LIBSUFFIX@.la ../mdlib/libmd@LIBSUFFIX@.la ../gmxlib/libgmx@LIBSUFFIX@.la

# AM_CPPFLAGS += -I/usr/local/openmm/include
# LDADD += /usr/local/openmm/lib/libOpenMM.*
# noinst_LTLIBRARIES += libopenmminterface.la
# libopenmminterface_la_SOURCES = md_openmm.cpp
# mdrun_LINK = $(CXXLINK)
# md1_LINK = $(CXXLINK)
# md2_LINK = $(CXXLINK)


libgmxpreprocess@LIBSUFFIX@_la_LIBADD       =
libgmxpreprocess@LIBSUFFIX@_la_DEPENDENCIES =

libgmxpreprocess@LIBSUFFIX@_la_SOURCES = \
add_par.c	add_par.h	\
compute_io.c	compute_io.h	\
convparm.c 	convparm.h	\
gen_ad.c 	gen_ad.h	\
gen_vsite.c	gen_vsite.h	\
genhydro.c	genhydro.h	\
gmxcpp.c	gmxcpp.h	\
gpp_atomtype.c		\
gpp_bond_atomtype.c	gpp_bond_atomtype.h	\
h_db.c		h_db.h		\
hackblock.c	hackblock.h	\
hizzie.c	hizzie.h	\
pdb2top.h	pdb2top.c	\
pgutil.c	pgutil.h	\
readir.c	readir.h	\
readpull.c	\
resall.c	resall.h	\
sorting.c 	sorting.h	\
specbond.c	specbond.h	\
ter_db.c	ter_db.h	\
tomorse.c	gpp_tomorse.h	\
topdirs.c	topdirs.h	\
topexcl.c	topexcl.h	\
topio.c 	topio.h		\
toppush.c 	toppush.h	\
topshake.c 	topshake.h	\
toputil.c	toputil.h	\
tpbcmp.c	tpbcmp.h	\
vsite_parm.c	vsite_parm.h	\
xlate.c		xlate.h

# There are several ways of adding independent programs that are
# independent of the package into Makefile.am
#
# METHOD 1: using autotools
# ===============================
#
# prog_SOURCES = prog.c header.h
# prog_LDADD =
#
# since prog has no dependency on others (except standard libraries),
# LDADD should be empty.  The advantage is that autotools will scan the
# source code and build the dependency for you.  However, sometimes this
# process can be unsuccessful, in which case .deps dir is empty, and it
# is unreliable in my opinion
#
# METHOD 2: direct compilation
# ================================
#
# prog$(EXEEXT): prog.c header.h
# 	$(COMPILE) -o $@ $<
#
# In this case, autotools will not try to generate dependencies for you,
# so it relies on the program make to do all the work.  This is the cheapest
# way, and I think is better.
#
# NOTE 0: $(EXEEXT) means .exe in windows, or nothing in linux/unix
#         $@ means whatever before the colon `:' on the first line,
#            i.e., prog$(EXEEXT) in this case
#         $< means the first file after the colon, i.e., prog.c
#
# NOTE 1: prog$(EXEEXT) is followed by a `:' for dependency, not an `=' as
# in automake.
#
# NOTE 2: make sure to use TAB, no space, at the beginning of the second line.
#
# NOTE 3: $(COMPILE) contains all exotic compiler flags that might make your
# program hard to compile or debug (e.g., INTEL C with MKL can lead to an
# unaddressable byte error with valgrind, while a plain compilation is free
# from this error).  Therefore, I believe sometimes its better to use $(CC)
# instead of $(COMPILE).  Several options that are included in $(COMPILE) are
#
# $(DEF) $(DEFAULT_INCLUDES) $(INCLUDES)
# $(AM_CPPFLAGS) $(CPPFLAGS)
# $(AM_CFLAGS) $(CFLAGS)
#
# Here, $(DEF) usually is just -DHAVE_CONFIG_H;  $(INCLUDES) is not the
# environment variable $INCLUDE, no `S' at the end; AM_  stands for automake.
#
# Another option is $(CC_FOR_BUILD), since some compilers only compile and
# do not produce an executable by default.
#
# It is desirable to create a macro SIMPLE_CC as

SIMPLE_CC = $(CC_FOR_BUILD) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
            $(AM_CPPFLAGS) $(CPPFLAGS) $(CFLAGS)

# Note: if there's a trouble remove $(CFLAGS) it may have exotic flags, e.g. -openmp
# Now for a simple program, we just use
#
# prog$(EXEEXT): prog.c header.h
# 	$(SIMPLE_CC) -o $@ $<

# On a side note, when you are using a VPATH (out-of-source) build, make sure
# that the source tree is clean.  Because if there is a prog$(EXEEXT) under
# the source directory, make will not even try to build a program in the build path.
# Therefore it is impossible to simultaneously use in-source and out-of-source
# build at the same time.

TOYS = ncontact dihmb1 dihmb2 amberize mkspx dihmb spbmb iibb \
       mqsub xdihs mkndx2 vdist mhex clusmb rmsdis rmsfilter mhidx

if MD_ONLY
bin_PROGRAMS = mdrun md1 md2 md2conv $(TOYS)
else
bin_PROGRAMS = \
	grompp		mdrun		tpbconv	  	pdb2gmx		\
	protonate	luck		gmxdump				\
	gmxcheck	x2top    \
	md1 md2 md2conv $(TOYS)
endif

mdrun_SOURCES = \
	glaasje.c 	glaasje.h 	gctio.c 	\
	ionize.c 	ionize.h 	xmdrun.h	\
	do_gct.c 	repl_ex.c	repl_ex.h	\
	xutils.c	md.c		mdrun.c		\
	genalg.c	genalg.h

# Theoretically, a single line like the following should work,
#
#   md1_SOURCES = md1.c
#
# There is a program depcomp that will take care of all the
# dependency issues.  However, but on some systems, depcomp
# is broken (e.g., biou), as one can see that the
# .deps directory is empty. this becomes problematic.
#
# The situation is especially trick because the inclusion
# is nested.
# md1.c <-- md1util.h <-- md1core.h <-- mb1.h, zcom1.h
#
# To resolve the program, we add a line to explicit declare
# the dependency of md1.$(OBJEXT) on these source files.
#
# Note: md1.$(OBJEXT) is followed by a colon `:' not an equal `='.
#
md1_SOURCES = md1.c
md1.$(OBJEXT): md1.c md1core.h md2util.h mb1.h zcom1.h

md2_SOURCES = md2.c
md2.$(OBJEXT): md2.c md2core.h md2util.h md2spb.h md2bb.h mb2.h zcom2.h

md2conv$(EXEEXT): md2conv.c md2spb.h zcom.h md2spbx.h
	$(SIMPLE_CC) -o $@ $< -lm

gmxdump_SOURCES = gmxdump.c

grompp_SOURCES = grompp.c

tpbconv_SOURCES = tpbconv.c

pdb2gmx_SOURCES = pdb2gmx.c

protonate_SOURCES = protonate.c

gmxcheck_SOURCES = gmxcheck.c

x2top_SOURCES = x2top.c nm2type.c x2top.h

mkndx2_SOURCES = mkndx2.c

dihmb_SOURCES = dihmb.c

spbmb_SOURCES = spbmb.c

#iibb_SOURCES = iibb.c
iibb$(EXEEXT): iibb.c
	$(SIMPLE_CC) -DGMXVERSION=40007 -o $@ $< -lm


ncontact_SOURCES = ncontact.c

dihmb1_SOURCES = dihmb1.c

dihmb2_SOURCES = dihmb2.c

vdist_SOURCES = vdist.c

mhex_SOURCES = mhex.c

clusmb_SOURCES = clusmb.c

rmsdis_SOURCES = rmsdis.c

rmsfilter_SOURCES = rmsfilter.c

mhidx_SOURCES = mhidx.c

# independent programs
amberize$(EXEEXT): amberize.c
	$(SIMPLE_CC) -o $@ $< -lm

mkspx$(EXEEXT): mkspx.c
	$(SIMPLE_CC) -o $@ $< -lm

mqsub$(EXEEXT): mqsub.c
	$(SIMPLE_CC) -o $@ $<

xdihs$(EXEEXT): xdihs.c zcom.h
	$(SIMPLE_CC) -o $@ $< -lm


# shortcut to install only mdrun

install-mdrun: mdrun$(EXEEXT)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	if test -f mdrun$(EXEEXT); then \
	  f=`echo mdrun$(EXEEXT)|sed 's/$(EXEEXT)$$//;$(transform);s/$$/$(EXEEXT)/'`; \
	  echo " $(INSTALL_PROGRAM_ENV) $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) mdrun$(EXEEXT) $(DESTDIR)$(bindir)/$$f"; \
	  $(INSTALL_PROGRAM_ENV) $(LIBTOOL)  --mode=install $(INSTALL_PROGRAM) mdrun$(EXEEXT) $(DESTDIR)$(bindir)/$$f; \
	else :; fi

CLEANFILES = *~ \\\#*

