#########################################
## Make file for arithmeticAverage-fsg ##
#########################################


##### Insert here location of headers and libraries
###################################################

# headers
ARMADILLO_INC	= /media/mircea/Stuff/Mircea/Programming/C++/libraries/armadillo-6.100.1/include

# libraries
OPENBLAS_LIB 	= /usr/lib/openblas/lib
#LAPACK_LIB		= /usr/lib/lapack

###################################################
###################################################

CC				= g++
DEBUG			= -DARMA_NO_DEBUG -DNDEBUG
PROG            = fsg
OBJECTS			= main.o binomialStrategy.o option.o
FILES			= 
INCLUDES        = -I $(ARMADILLO_INC)
LIBS            = -L $(OPENBLAS_LIB) -lopenblas
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
arithmeticAverage-fsg.o: helpers.h fsg_average.h binomialStrategy.h option.h
#
binomialStrategy.o: binomialStrategy.h
#
option.o: option.h


clean:
	\rm *.o *~ fsg
