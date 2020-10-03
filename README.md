# StatisticalMechanics
#Los diferentes códigos aquí alojados son escritos en lenguaje c++, empleando las librerías de root/cern, y gsl,  para compilar y ejecutar el fichero Ising.cpp, se emplea el comando en terminal

g++ Ising.cpp `root-config --glibs --cflags` -c -o asd.o; g++ asd.o `root-config --glibs --cflags` -o Ising.exec; ./Ising.exec
