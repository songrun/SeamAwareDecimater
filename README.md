# Seamless

## Interesting article: http://rastergrid.com/blog/2010/09/history-of-hardware-tessellation/
## Compiling

This code depends on:

- [libigl](https://github.com/libigl/libigl)
- [eigen](http://eigen.tuxfamily.org/) (e.g. `brew install eigen`)
- [GLFW3](http://www.glfw.org/) (e.g. `brew install glfw3`)

### Download libigl and compile the third-party dependencies
    git clone https://github.com/libigl/libigl.git --recursive

### Compile this project
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make
