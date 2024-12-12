# C-GLASS notes





## C-GLASS flowchart

## TriMesh
### Initialization flowchart

From `ply` file:
```
Init(system_parameters *params)
├── SetParameters()
│   ├── rng_
│   ├── ply_path
│   └── ...phsyical params...
├── LoadPly()
│   └── ...he matrices...
├── SyncWithMats()
│   ├── vrts_::index_,pos_,
│   ├── edges_
│   ├── tris_
│   ├── half_edges_
│   └── boundaries_
├── DrawVerts()
├── UpdateIncidenceData()
│   ├── v.UpdateNeighbors() for v in vrts_
│   │   ├── neighbs_, edges_, tris_
│   │   └── n_neighbs_, n_edges_, n_tris_
│   ├── e.UpdateIncidentTris() for e in edges_
│   └── f.UpdateIncidentEdges() for f in tris_
└── UpdateGeometricData()
    ├── e.UpdateEdgeGeometry() for e in edges_
    │   └──  vector_, length_
    ├── f.UpdateArea(), f.UpdateNormal() for f in tris_
    ├── UpdateCentroid()
    └── f.UpdateVolume(centroid_) for f in tris_
```

From icosohedron:
```
Init(system_parameters *params)
├── SetParameters()
│   ├── rng_
│   ├── ply_path
│   └── ...phsyical params...
├── MakeIcosphere()
│   └── ...
├── InitializeMeshOG()
│   ├── vrts_::index_,pos_,
│   ├── edges_
│   ├── tris_
│   ├── half_edges_
│   └── boundaries_
|
└──...
```

### Time evolution flowchart


## Table of Contents
- [Installation](#installation)
  - [Requirements](#requirements)
    - [Arch](#arch)
    - [Equivalent Ubuntu/Arch packages](#equivalent-ubuntuarch-packages)
  - [Using C-GLASS with Virtualbox](#using-c-glass-with-virtualbox)
    - [EndeavourOS](#endeavouros)
  - [Install script](#install-script)
    - [Compilation error in logger.cpp](#compilation-error-in-loggercpp)
    - [Help CMake find glew](#help-cmake-find-glew)
- [Running C-GLASS](#running-c-glass)
- [To do](#to-do)

## Installation

Clone with submodules
```bash
git clone -b nab_xfer --recursive https://github.com/wlough/C-GLASS
```

### Requirements

#### Arch
Most required packages (I'm looking at you `armadillo`) are in the official Arch repos:
```bash
sudo pacman -S \
cmake \
gsl \
openmpi \
fftw \
yaml-cpp \
boost-libs \
gcc \
pkgconf \
glfw \
glew \
boost
```
maybe need:
```bash
sudo pacman -S \
libxi \
libxcursor \
libxinerama \
```
from aur:
```bash
yay -S \
armadillo
```

<!-- 
```bash
sudo pacman -S \
git \
npm \
clang-tools-extra \
cmake \
gsl \
openmpi \
fftw \
yaml-cpp \
boost-libs \
gcc \
pkgconf \
glfw-x11 \
glew \
libxi \
libxcursor \
libxinerama \
boost
``` -->

#### Equivalent Ubuntu/Arch packages
```
cmake
libyaml-cpp-dev --> yaml-cpp
libgsl-dev --> gsl
libopenmpi-dev --> openmpi
libfftw3-dev --> fftw
libboost-math-dev --> boost-libs
g++ --> gcc
libarmadillo-dev --> armadillo
pkg-config --> pkgconf
libglfw3-dev --> glfw-x11
libglew-dev --> glew
libxi-dev --> libxi
libxcursor-dev --> libxcursor
libxinerama-dev --> libxinerama
```

### Using C-GLASS with Virtualbox

- [How to run an Ubuntu virtual machine](https://ubuntu.com/tutorials/how-to-run-ubuntu-desktop-on-a-virtual-machine-using-virtualbox)
- [How to install Virtualbox on EndeavourOS](https://discovery.endeavouros.com/applications/how-to-install-virtualbox/2021/03/)
- [How to run an EndeavourOS virtual machine](https://discovery.endeavouros.com/applications/installing-endeavouros-on-virtualbox/2021/03/)

#### EndeavourOS

Add shared folder to guest machine accesible to user `vuser`:
- make shared directory on host machine: `mkdir /host/path/to/share`
- You don't need to explicitly make directory on guest machine
- Add shared folder in virtualbox settings
  - Folder path = `/host/path/to/share`
  - Mount point = `/guest/path/to/share`
- On guest machine, add `vuser` to `vboxsf` group with `sudo usermod -aG vboxsf vuser`
  - members of `vboxsf` group should already have read/write privileges for `share`, but you can check with `sudo ls -l /guest/path/to/share`

```bash
mkdir ~/share
sudo mount -t vboxsf share ~/share/
```


### Install script

```bash
./install.sh -c && ./install.sh -wg
```
With graphics:

```bash
./install.sh -wg
```
Clear previous install files
```bash
./install.sh -c
```

#### Compilation error in logger.cpp

If you see the error

```bash
[ 49%] Building CXX object src/CMakeFiles/cglass.dir/logger.cpp.o
/home/wlough/git/C-GLASS/src/logger.cpp: In static member function ‘static void Logger::Error(const char*, ...)’:
/home/wlough/git/C-GLASS/src/logger.cpp:128:3: error: ‘exit’ was not declared in this scope
  128 |   exit(1);
      |   ^~~~
/home/wlough/git/C-GLASS/src/logger.cpp:2:1: note: ‘exit’ is defined in header ‘<cstdlib>’; did you forget to ‘#include <cstdlib>’?
    1 | #include "cglass/logger.hpp"

```

add the `#include <cstdlib>` directive at the top of the `C-GLASS/src/logger.cpp`:

```C
#include <cstdlib>  // add this line

#include "cglass/logger.hpp"

/****************************/
/******** SINGLETON *********/
/****************************/
...
```

#### Help CMake find glew

If you see the error

```bash
CMake Error at src/CMakeLists.txt:23 (find_package):
  By not providing "Findglew.cmake" in CMAKE_MODULE_PATH this project has
  asked CMake to find a package configuration file provided by "glew", but
  CMake did not find one.

  Could not find a package configuration file provided by "glew" with any of
  the following names:

    glewConfig.cmake
    glew-config.cmake

  Add the installation prefix of "glew" to CMAKE_PREFIX_PATH or set
  "glew_DIR" to a directory containing one of the above files.  If "glew"
  provides a separate development package or SDK, be sure it has been
  installed.
```

then try renaming

`C-GLASS/.CMake_Modules/FindGLEW.cmake`

file to

`C-GLASS/.CMake_Modules/Findglew.cmake`







## Running C-GLASS
```bash
cd C-GLASS
./bin/build/executable/cglass.exe params.yaml
```









## To do
* Define primitive interactor-elements `meshbrane::primitives`
  * `Point`, `Vertex`, `OrientedSegment`, `OrientedTri`,...?
* Define composite objects
  * OrientedSurface
    * RigidBody.OrientedSurface
    * Membrane(OrientedSurface)
* Define species parameters
  * Add to `config/default_config.yaml` and run `./configure_cglass.exe config/default_config.yaml`
* Define interaction potentials: tethering potential
  * Define subclass of PotentialBase
  * `#include NewPotential` in PotentialManager
  * Add new potential_type to definitions.hpp for lookup purposes
  * See InitPotentials method in PotentialManager.h for examples
<!-- * Define output types for new species (posit, spec, checkpoint) -->
