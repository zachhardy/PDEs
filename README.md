# PDEs
A C++ framework for numerically solving PDEs.

## Compiling
 
The **PDEs** library requires an MPI library (**MPICH** or **OpenMPI**) and the 
linear algebra package **PETSc**. Packages such as **HomeBrew** and **MacPorts** may be
used to install these, but for manual installation, these packages can be found at:
<ul>
    <li> <b>MPICH</b>: https://www.mpich.org/ </li>
    <li> <b>OpenMPI</b>: https://www.open-mpi.org/ </li>
    <li> <b>PETSc</b>: https://petsc.org/release/ </li>
</ul>

### Installing MPI

From the **MPICH** or **OpenMPI** site, navigate to the releases, download an appropriate 
version, and extract it. Once extracted, set up a build and install directory. Normally 
the install directory is located outside the source. Navigate to the build directory, 
then call:

```console
    $ /path/to/source/configure -prefix=/path/to/mpi-install
    $ make
    $ make install
```

For more options or troubleshooting, refer to the documentation of the MPI 
package. 

Once this is completed, the path to the executables should be added to the
path via

```command
    $ export PATH=/path/to/mpi-install/bin:$PATH 
```

Note that if this is not within the terminal startup script 
(<tt>.bash_profile</tt>, <tt>.bashrc</tt>, <tt>.zshrc</tt>, etc.) this will 
have to be called each time a new terminal session is started.

If the installation was successful, the command 

```command 
    $ mpicc --version
``` 

should yield  a result which displays the type and version of the MPI compiler. 
If a different  MPI version is referenced than expected there are a number of
options. First, one can use an alias:

```command
    $ alias mpicc=/path/to/mpi-install/bin/mpicc
```

Alternatively,  the terminal startup script can be modified to ensure that the 
path to the correct MPI compilers appear before any other. Lastly, one can 
always use the full path to the compiler to use it.

If issues persist, consult the documentation from the MPI documentation.

### Installing PETSc

From the **PETSc** website, navigate to the relases, download an appropriate version,
and extract it. Alternatively, one can use ```wget``` via:

```bash
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-<version>.tar.gz
tar -zxf petsc-<version>.tar.gz
```

Once extracted, create an install directory and navigate to the source directory
<tt>petsc-&lt;version&gt;</tt>. PETSc should be configured using the following options:

```bash
./configure  \
--prefix=/path/to/petsc-install \
--download-hypre=1  \
--with-ssl=0  \
--with-debugging=0  \
--with-pic=1  \
--with-shared-libraries=1  \
--download-fblaslapack=1  \
--download-metis=1  \
--download-parmetis=1  \
--download-superlu_dist=1  \
--with-cxx-dialect=C++11  \
CFLAGS='-fPIC -fopenmp'  \
CXXFLAGS='-fPIC -fopenmp'  \
FFLAGS='-fPIC -fopenmp'  \
FCFLAGS='-fPIC -fopenmp'  \
F90FLAGS='-fPIC -fopenmp'  \
F77FLAGS='-fPIC -fopenmp'  \
COPTFLAGS='-O3 -march=native -mtune=native'  \
CXXOPTFLAGS='-O3 -march=native -mtune=native'  \
FOPTFLAGS='-O3 -march=native -mtune=native'  \
PETSC_DIR=$PWD
```

If the configuration fails, consult the **PETSc** documentation. If 
successful, execute the <tt>make</tt> command provided by PETSc.  After 
**PETSc** finishes the build, execute the provided <tt>make install</tt> 
command. Once this is completed, the install can be tested using
```make test```. 

As before, an environment variable for the PETSc install directory should 
be set via:

```command
    $ export PETSC_ROOT=/path/to/petsc-install
```

As before, it is recommended that this be included in the terminal startup
script.


### Documentation

The documentation can be generated locally via:

```command
    $ doxygen doc/DoxyConfig
```

and can be viewed via HTML using the file <tt>doc/html/index.html</tt>.
