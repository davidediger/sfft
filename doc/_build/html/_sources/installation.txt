Installation
============

Prerequisites
-------------

The SFFT library was only tested on Linux systems and is only guaranteed to
work there. However, the library should also be able to compile on other 
platforms and operating systems.

The following packages have to be installed to compile the library:

- Python (any version > 2.3, including Python 3), used by the waf_ build
  system
- FFTW_ 3
- *(optionally)* `Intel Integrated Performance Primitives`_

If you want to build benchmark tools, also install

- Valgrind_

The SFFT library is known to work the following compilers:

- GCC_ (tested with GCC 4.4 and 4.7)
- `Intel C++ Compiler`_  (only versions >= 13, does NOT work with ICC 12)



.. _waf: https://code.google.com/p/waf/
.. _FFTW: http://www.fftw.org/
.. _Intel Integrated Performance Primitives: 
         http://software.intel.com/intel-ipp
.. _Valgrind: http://valgrind.org/ 
.. _GCC: http://gcc.gnu.org/
.. _Intel C++ Compiler: http://software.intel.com/en-us/intel-compilers



Compiling From Source and Installation
--------------------------------------

Unpack the tarball and change into the newly created directory 
(sfft-*version*). Then, the SFFT library can be built with a simple::

    $ ./configure
    $ make

and installed with::

    $ make install

Some configuration options can be passed to the configuration script. The most
important are::

    $ ./configure --help

    [...]
    --debug               compile in debug mode
    --profile             add source-level profiling to instruction counting programs
    --without-ipp         do not the Intel Performance Primitives library
    [...]

Use ``--debug`` and ``--profile`` are only useful when developing (see
:ref:`development`). The option ``--without-ipp`` is to be used when you do not have
Intel IPP installed.

When these steps succeded, you should be ready to use the SFFT library.


Linking against the SFFT Library
--------------------------------

Two versions of the SFFT library are built when compiling the sourcecode:
a static library (libsfft.a) and a shared library (libsfft.so). You can link
these libraries in your programs like any other library, but you have to make
sure that you link dependencies as well.

Do not forget to link:

- FFTW, for example via *pkg-config*: ``pkg-config --cflags --libs fftw3``
- Intel IPP (if not disabled via ``--without-ipp``), 
  e.g. ``-lippvm -lipps -pthread``
- Your compilers OpenMP library, for example ``-lgomp`` for GCC
- *libm* and *librt* (``-lm -lrt``)
