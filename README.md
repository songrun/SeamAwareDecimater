# Seam-aware Decimater

This project implements the seam-aware decimation portion of our SIGGRAPH Asia 2017 paper 
[Seamless: Seam erasure and seam-aware decoupling of shape from mesh resolution](https://cragl.cs.gmu.edu/seamless/).
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
