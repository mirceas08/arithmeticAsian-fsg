###############################################
## Make file for Heston model discretization ##
###############################################


##### Insert here location of headers and libraries
###################################################

# headers
ARMADILLO_INC	= /media/mircea/Stuff/Mircea/Programming/C++/libraries/armadillo-6.100.1/include

# libraries
OPENBLAS_LIB 	= /usr/lib/openblas/lib
LAPACK_LIB		= /usr/lib/lapack

###################################################
###################################################

CC				= g++
DEBUG			= -DARMA_NO_DEBUG -DNDEBUG
PROG            = arithmeticAverage-fsg
OBJECTS			= arithmeticAverage-fsg.o
FILES			= 
INCLUDES        = -I $(ARMADILLO_INC)
LIBS            = -L $(OPENBLAS_LIB) -lopenblas -L $(LAPACK_LIB) -llapack
CFLAGS        	= -O3 -std=c++11 -fopenmp -march=native $(DEBUG)

arithmeticAverage-fsg: $(OBJECTS)
	$(CC) -o $(PROG) $(OBJECTS) $(INCLUDES) $(LIBS) $(CFLAGS)


## Generic rules:
.SUFFIXES:          # get rid of that annoying Modula rule
.SUFFIXES: .cpp .h .o

#
# Rules
#
.cpp.o:
	$(CC) -c $*.cpp $(INCLUDES) $(LIBS) $(CFLAGS)
#
# Dependencies
#
arithmeticAverage-fsg.o: helpers.h forward_shooting_grid.h fsg_vector.h


clean:
	\rm *.o *~ arithmeticAverage-fsg
