#NAME=PC
NAME=LAPTOP
ifeq ($(NAME), PC)
INC=-I/usr/local/Cellar/boost/1.58.0/include -I/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/include -I/usr/local/include -I/Users/juncheng/work/phosim_core/source/cfitsio/include -I/usr/include/python2.7/
LIB=-L/usr/local/Cellar/boost/1.58.0/lib -L/usr/local/lib -L/Users/cheng109/toberemoved/phosim/phosim_core/source/cfitsio/lib -L/opt/local/lib
INC_EIGEN=-I/Users/cheng109/eigen-eigen-c58038c56923/
endif
ifeq ($(NAME), LAPTOP)
INC=-I/usr/local/Cellar/boost/1.58.0/include -I/Users/juncheng/work/cfitsio -I/usr/local/include -I/usr/local/Cellar/gsl/1.16/include
LIB=-L/usr/local/Cellar/boost/1.58.0/lib -L/Users/juncheng/work/cfitsio -L/opt/local/lib -L/usr/local/Cellar/gsl/1.16/lib
INC_EIGEN=-I/Users/juncheng/work/eigen-eigen-b30b87236a1b
endif
CC=g++ #clang++
CFLAGS=-Wall -O3 $(INC) $(INC_EIGEN) $(LIB)
LDFLAGS=-lcfitsio -lgsl -L. -lfortranstuff 
#-lgfortran 
#-lg2c
# -larmadillo -lboost_iostreams -lboost_system #-fopenmp 
SRCS=main.cpp Image.cpp commons.cpp Model.cpp
FC=gfortran

COMMON_HDRS= Image.h Model.h commons.h fastell.h gl_crit.h parafit.h  
OBJS= commons.o junGL.o Image.o Model.o gl_crit.o parafit.o
#OMP_NUM_THREADS=8

#all: main.o Image.o commons.o Model.o
#	g++ $(INC) $(LIB) $(CFLAG) $(LDFLAGS) commons.cpp main.cpp Image.cpp Model.cpp -o junGL 
#	./junGL
all: $(COMMON_HDRS) fastell_example.o

	$(CC) $(CFLAGS) -o junGL $(LDFLAGS) commons.cpp main.cpp Image.cpp Model.cpp gl_crit.cpp parafit.cpp #-mmacosx-version-min=10.8
	#valgrind --tool=memcheck --leak-check=full --verbose --log-file=memcheck.log --track-origins=yes ./junGL
	./junGL

#libfortranstuff.a: 
#	$(FC) -O -c slatec/src/*.f fastell.f
#	ar -r libfortranstuff.a *.o
#	rm *.o

fastell_example.o: 
	$(FC) -c fastell_example.f 
	$(CC) -c testFastEll.cpp 

plot:
	#@export column=5
	#@export fileName=output.txt
	./plotScript
#fastell.o: 
#Model.o: 
#	$(CC) $(CFLAGS) $(LDFLAGS) -c Model.cpp 


#test:
#	$(FC) slatec/src/*.f fastell_example.f -o testFastEll
#	./testFastEll