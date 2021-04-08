#
# Makefile for CBS program
# Emanuele Bosoni 10/01/2018
#
# The idea is to use the fdf directory in the Siesta Scr and his arch.make
# This is achieved by the VPATH directive below.
# The arch.make file is supposed to be in $(OBJDIR). This is normally
# the top Obj, but if you are using architecture-dependent build directories
# you might want to change this. (If you do not understand this, you do not
# need to change anything. Power users can do "make OBJDIR=Whatever".)
#
# As the program is serial, be sure to define an FC_SERIAL symbol in 
# your top arch.make. Lapack/Blas routines are nedded, the ones 
# in arch.make are used as default
#
OBJDIR=Obj
#
.SUFFIXES: .f .F .o .a .f90 .F90
#
VPATH:=$(shell pwd)/../../Src
#
default: what
#
ARCH_MAKE=../../$(OBJDIR)/arch.make
include $(ARCH_MAKE)
#
FC_DEFAULT:=$(FC)
FC_SERIAL?=$(FC_DEFAULT)
FC:=$(FC_SERIAL)         # Make it non-recursive
#
DEFS:=$(DEFS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
FPPFLAGS:=$(FPPFLAGS) $(DEFS_PREFIX)-UMPI $(DEFS_PREFIX)-UCDF
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#
SCALAPACK=$(SCALAPACK_LIBS) $(LAPACK_LIBS) $(BLACS_LIBS) $(BLAS_LIBS)

what:
	@echo "Using the arch.make in $(ARCH_MAKE)"
	@echo "Compilation architecture to be used: ${SIESTA_ARCH}"
	@echo "If this is not what you want, create the right"
	@echo "arch.make file using the models in Sys"
	@echo
	@echo "Type make cbs to compile"
#@sleep 2
#
SYSOBJ=$(SYS).o
#
# Note that machine-specific files are now in top Src directory.
#
OBJS_STM = hsx_m.o calccbs.o cbs.o

INCFLAGS:= $(NETCDF_INCFLAGS) $(INCFLAGS)
#
# Use the makefile in Src/fdf and all the sources there.
#
FDF=libfdf.a
FDF_MAKEFILE=$(VPATH)/fdf/makefile
FDF_INCFLAGS:=-I $(VPATH)/fdf $(INCFLAGS)
$(FDF): 
	(mkdir -p fdf ; cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "FC=$(FC)" "VPATH=$(VPATH)/fdf" \
                          "DEFS=$(DEFS)" \
                          "FPPFLAGS=$(FPPFLAGS)"\
                          "ARCH_MAKE=../$(ARCH_MAKE)" \
                          "INCFLAGS=$(FDF_INCFLAGS)" "FFLAGS=$(FFLAGS)" module)
#

COM_OBJS_STM=$(OBJS_STM) $(SYSOBJ)
ALL_OBJS_STM=$(MOD_OBJS) $(COM_OBJS_STM)
#
cbs: $(FDF)  $(ALL_OBJS_STM)
	$(FC) -o cbs \
	       $(LDFLAGS) $(ALL_OBJS_STM) $(FDF) $(SCALAPACK) 
#
clean:
	@echo "==> Cleaning object, library, and executable files"
	rm -f cbs *.o  *.a *.mod
	(cd fdf ; $(MAKE) -f $(FDF_MAKEFILE) "ARCH_MAKE=../$(ARCH_MAKE)" clean)
	rm -f _tmp_deps deps.list  protomake*
#
#dep:
#	rm -f _tmp_deps deps.list  protomake*
#	sfmakedepend --depend=obj --file _tmp_defs  --modext=o \
#          $(VPATH)/*.f $(VPATH)/*.f90 $(VPATH)/*.F $(VPATH)/*.F90 \
          *.f *.F *.f90 *.F90
#	@echo "Ignore errors in the following command"
#	@echo "They appear if the last grep does not match anything"
#	@-for i in $(ALL_OBJS_STM) ; do grep "^$$i: " _tmp_defs  ; done > deps.list
#	@echo "ok"
#	@csplit -k -f "protomake" Makefile '/^# -- specially prepared/+1'
#	@cat protomake00 deps.list >| Makefile
