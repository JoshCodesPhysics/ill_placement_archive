#SHELL = /bin/posix/sh 
# derni�re version qui marche 12 ou 12.1 
# Version 12
#SRCF = lett_min_x.f90 \
       lect_ligne_x.f90 \
       lect_entree_x.f90 \
       invmat.f90 \
       fctmu_1.f90 \
       env_12.f90 \
       coordrot_2.f90 
# Version 12.1
#SRCF = lett_min_x.f90 \
       lect_ligne_x.f90 \
       lect_entree_x.f90 \
       invmat.f90 \
       fctmu_1.f90 \
       env_12.1.f90 \
       coordrot_2.f90 


# sources
# Version 14
#SRCF = lett_min_x.f90 \
#       lect_ligne_x.f90 \
#       lect_entree_x.f90 \
#       case.f90 \
#       invmat.f90 \
#       fctmu_1.f90 \
#       env_14.f90 \
#       coordrot_2.f90 
#
# sources
# Version 15 
SRCB = lett_min_x.f90 \
       lect_ligne_x.f90 \
       lect_entree_x.f90 \
       case.f90 

SRCF = sym3.f90	\
       invmat.f90 \
       fctmu_1.f90 \
       coordrot_2.f90  \
       env_15.f90 \
       sym.f90

SRCT = test_sym.f90 


## Version 15 symetrie D2h testee

# objets
OBJF = $(SRCF:.f90=.o)
OBJB = $(SRCB:.f90=.o)
OBJT = $(SRCT:.f90=.o)

# executables
PROG = xenv15
XTEST = xtest


#------------------------------------------------------
#  Ordinateur de Marie + calc-maitre
# Compilateurs
F90 = gfortran     #ifort
CC  = gcc          #icc
FFLAGS= -O2 
#FFLAGS= -g  -fbounds-check
CFLAGS =    
LFLAGS=       # -Bstatic
LDLIBS =      # -Bstatic
#------------------------------------------------------
#------------------------------------------------------
#  Crihan + idriss 
# Compilateurs
#F90 = ifort
#CC  = icc
#FFLAGS= -O2 
#FFLAGS= -g  -fbounds-check
#CFLAGS =    
#LFLAGS=        -Bstatic
#LDLIBS =       -Bstatic
#------------------------------------------------------

# drapeaux
#LFLAGS = -Bstatic -\# -C -O3 -R/usr/lib64/gcc/x86_64-suse-linux/4.1.2
#FFLAGS = -fast -tp px



main :  $(OBJB) $(OBJF)
	$(F90) $(LFLAGS) $(OBJF) $(OBJB) $(LDLIBS) -o $(PROG) 

test :  $(OBJB) $(OBJT) 
	$(F90) $(LFLAGS) $(OBJT) $(OBJB) $(LDLIBS) -o $(XTEST) 

$(OBJF): $(SRCF) 
	 $(F90) -c $(FFLAGS) $(SRCF) 

$(OBJB): $(SRCB) 
	 $(F90) -c $(FFLAGS) $(SRCB) 

$(OBJT): $(SRCT) 
	 $(F90) -c $(FFLAGS) $(SRCT) 

clean:
	@echo Removing object files ... 
	rm -f $(OBJC) $(OBJF) $(OBJB) $(OBJT) *.mod



