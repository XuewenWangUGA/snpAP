
all: Bam2snpAD bamtools snpADJoin

bamtools/build/src/libbamtools.a: 
	rm -rf bamtools
	git clone https://github.com/pezmaster31/bamtools.git
	mkdir bamtools/build || true
	cd bamtools/build ; cmake .. ; make

Bam2snpAD: bamtools/build/src/libbamtools.a src/BamAliIt.cpp src/BamTable.cpp src/FastaSlurp.cpp 
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -I bamtools/build/src -I bamtools/src -c -o src/BamAliIt.o src/BamAliIt.cpp 
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -I bamtools/build/src -I bamtools/src -c -o src/BamTable.o src/BamTable.cpp 
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -I bamtools/build/src -I bamtools/src -c -o src/FastaSlurp.o src/FastaSlurp.cpp 
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -I bamtools/build/src -I bamtools/src -L bamtools/build/src -o Bam2snpAD src/BamAliIt.o src/BamTable.o src/FastaSlurp.o -lpopt -lz -lbamtools -lstdc++ -lz

snpADJoin: src/snpADJoin.cc src/SnpADFile.cc src/SnpADFile.h
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -c -o src/SnpADFile.o src/SnpADFile.cc
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -c -o src/snpADJoin.o src/snpADJoin.cc
	$(CXX) -std=c++11 -O3 -g -I src/ -Wall -o snpADjoin src/SnpADFile.o src/snpADJoin.o
	

clean:
	rm -rf src/*.o Bam2snpAD

super-clean:
	rm -f bamtools/lib/libbamtools.a
