# Seamless Decimater

This code depends on:

- [libigl](https://github.com/libigl/libigl)
- [eigen](http://eigen.tuxfamily.org/) (e.g. `brew install eigen`)

### Download libigl and compile the third-party dependencies
    git clone https://github.com/libigl/libigl.git --recursive

### Compile this project
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make
    
### Run this project
	./decimater ../models/animal.obj percent-vertices 50
	./decimater ../models/animal.obj num-vertices 1000
