CC := mpiCC
CPPFLAGS := -O3 -Wno-unused-result -Wno-format

amrmake: main.cpp
	$(CC) -o mpi_mesh main.cpp $(CPPFLAGS)
