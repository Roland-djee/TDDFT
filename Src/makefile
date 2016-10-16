#############################################################################
#                                                                           #
#	MAKEFILE GENERIQUE POUR LA COMPILATION DES CODES FORTRAN            #
#                                                                           #
#############################################################################

# Indiquer quel compilateur est à utiliser

COMP = ifort
#COMP = f90
#COMP = gfortran

# Spécifier les options du compilateur CFLAG = Compilation FLAGS,
# WFLAGS = Warning FLAGS, LDFLAGS = Link eDition FLAGS, OPTFLAGS = OPTimizing
# FLAGS, DEBUGFLAGS = DEBUGging FLAGS, LAPACK = édition de liens avec la 
# librairie LAPACKS, PROFLAGS = PROfiling FLAGS, PARFLAGS = PARallel FLAGS

# ifort FLAGS

#WFLAGS  = -warn argument checking -check underflow
#WFLAGS  = -C -check all -warn all,nodec,interfaces -gen_interfaces -traceback # available with ifort 
#WFLAGS  = -check all -warn all -traceback
#WFLAGS  = +check=all
#CFLAGS  = -c -g -heap-arrays # available with ifort
LAPACK = -L/usr/lib -llapack -lblas -lgfortran
#LAPACK = -llapack
#BLAS = -L -s /opt/fortran90/share/man/man3.Z/ -lblas
#PARFLAGS = +Oopenmp
PARFLAGS = -openmp
#OPTFLAGS = -O2 -mtune=core2
OPTFLAGS = -O3
#DEBUGFLAGS = -g # available with ifort

# gfortran FLAGS

#LAPACK = -llapack -lblas
#LAPACK = -llapack 
#PARFLAGS = -fopenmp

# Common FLAGS

#CFLAGS  = -c
CFLAGS = -c -g 
LDFLAGS = -o
#PROFLAGS = -pg

# Reconnaître les Extension de nom de fichier *.out, *.o, .f90 et .h
# comme suffixe

SUFFIXES = .exe .o .f90 
.SUFFIXES: $(SUFFIXES) .

# Nom de l'exécutable

EXEC = Atlas

# Liste des fichiers modules nécessaires pour le programme final

MOD = IonProbaLib.f90\
	Expokit.f90\
	Mathlib.f90\
	Library.f90

# Listes de fichiers sources nécessaires pour le programme final

SRC  = Atlas.f90 

# Liste de fichiers objets nécessaires pour le programme final.
# Il s'agit des mêmes objets que ceux source avec une extension
# .o.

OBJS  = $(MOD:.f90=.o)\
	$(SRC:.f90=.o)

all: $(EXEC) Class Mrproper
 
# Etape de compilation de tous les fichiers modules .f90 en .o

Modules:
	@$(COMP) $(CFLAGS) $(PARFLAGS) $(WFLAGS) $(MOD) $(PROFLAGS)

# Etape de compilation de tous les fichiers sources .f90 en .o

Sources:
	@$(COMP) $(PARFLAGS) $(CFLAGS) $(WFLAGS) $(OPTFLAGS) $(SRC) $(PROFLAGS)

# Etape d'édition de liens

Atlas: Modules Sources 
	@$(COMP) $(DEBUGFLAGS) $(OBJS) $(LAPACK) $(BLAS) $(PARFLAGS) $(LDFLAGS) $(EXEC).exe $(PROFLAGS)
	@echo "----Compilation done----"

# Classement des fichiers

Class:
	@mv *.o ../Objets/.
	@mv *.exe ../Exec/.
#	@mv *.f90~ ../Save/.
	@echo "----Data classification done----"

# Nettoyage

Clean:
	@rm -rf *.mod

Mrproper: Clean
	@rm -rf $(EXEC)
	@echo "----Data cleaning done----"

