OBJS = src/matrix-elements-dq.f90

FC = gfortran

COMPILER_FLAGS = -ffree-line-length-none -Wall -Wextra -g -ffloat-store -std=legacy -O3 -no-pie
LINKER_FLAGS = lib/drakelib-dq.a lib/dqlib.a -J modules
OBJ_NAME = build/matrix-elements-dq

all: 
	$(FC) $(OBJS) $(COMPILER_FLAGS) $(INCLUDE_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)

