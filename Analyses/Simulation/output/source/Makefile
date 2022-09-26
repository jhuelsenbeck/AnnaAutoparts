LIBS     =   -lm

CFLAGS   =  -fast 

CC       =   g++

OBJECTS  =   main.o Alignment.o Chunk.o Hungarian.o MbEigensystem.o MbMath.o MbRandom.o MbTransitionMatrix.o Mcmc.o Model.o Parm.o ParmAsrv.o ParmFreqs.o ParmLength.o ParmSubrates.o ParmTree.o Partition.o Restaurant.o Settings.o StateSets.o Table.o

PROGS    = AutoParts

AutoParts:		$(PROGS)

AutoParts:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o AutoParts
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Alignment.o:	Alignment.cpp
		$(CC) $(CFLAGS) -c Alignment.cpp

Chunk.o:	Chunk.cpp
		$(CC) $(CFLAGS) -c Chunk.cpp

Hungarian.o:	Hungarian.cpp
		$(CC) $(CFLAGS) -c Hungarian.cpp

MbEigensystem.o:	MbEigensystem.cpp
		$(CC) $(CFLAGS) -c MbEigensystem.cpp

MbMath.o:	MbMath.cpp
		$(CC) $(CFLAGS) -c MbMath.cpp

MbRandom.o:	MbRandom.cpp
		$(CC) $(CFLAGS) -c MbRandom.cpp

MbTransitionMatrix.o:	MbTransitionMatrix.cpp
		$(CC) $(CFLAGS) -c MbTransitionMatrix.cpp

Mcmc.o:	Mcmc.cpp
		$(CC) $(CFLAGS) -c Mcmc.cpp

Model.o:	Model.cpp
		$(CC) $(CFLAGS) -c Model.cpp

Parm.o:	Parm.cpp
		$(CC) $(CFLAGS) -c Parm.cpp

ParmAsrv.o:	ParmAsrv.cpp
		$(CC) $(CFLAGS) -c ParmAsrv.cpp

ParmFreqs.o:	ParmFreqs.cpp
		$(CC) $(CFLAGS) -c ParmFreqs.cpp

ParmLength.o:	ParmLength.cpp
		$(CC) $(CFLAGS) -c ParmLength.cpp

ParmSubrates.o:	ParmSubrates.cpp
		$(CC) $(CFLAGS) -c ParmSubrates.cpp

ParmTree.o:	ParmTree.cpp
		$(CC) $(CFLAGS) -c ParmTree.cpp

Partition.o:	Partition.cpp
		$(CC) $(CFLAGS) -c Partition.cpp

Restaurant.o:	Restaurant.cpp
		$(CC) $(CFLAGS) -c Restaurant.cpp

Settings.o:	Settings.cpp
		$(CC) $(CFLAGS) -c Settings.cpp

StateSets.o:	StateSets.cpp
		$(CC) $(CFLAGS) -c StateSets.cpp

Table.o:	Table.cpp
		$(CC) $(CFLAGS) -c Table.cpp

clean:		
		rm -f *.o
