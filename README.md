# Seam-aware Decimater

Seam-aware decimater simplifies a mesh while preserving its UV parameterization. It allows 
the same texture to be used across all decimation levels—notably along seams.

### Requirements

This project uses C++ 11, and it depends on:

- [libigl](https://github.com/libigl/libigl) (`git clone https://github.com/libigl/libigl.git --recursive`)
- [eigen](http://eigen.tuxfamily.org/) (e.g. `brew install eigen`)

### Compile this project
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make
    
### Run this project
	./decimater ../models/animal.obj percent-vertices 50
	./decimater ../models/animal.obj num-vertices 1000
