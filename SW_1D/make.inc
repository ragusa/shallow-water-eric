ROOT  = /users/eric/shallow-water-eric/SFeMANS
OPT   = -O -Wall -fimplicit-none -fbounds-check -g -C
#OPT   = -march=native -mtune=native -Ofast
F90   = gfortran $(OPT)
LNK   = gfortran $(OPT)

RACINE_INTEL = /w/guermond
LIB_PARDISO =
INC_PARDISO =

SPLIB_DIR =
SPLIB_OBJ =
