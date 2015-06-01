Usage
=====

All types and functions of the SFFT library are defined in the header
``sfft.h``. Include it at the beginning of your program.


Computing Sparse DFTs
---------------------

Creating Plans
++++++++++++++

SFFT executions consist of two seperate steps: planning and execution. The
planning phase is only executed once for specific input parameters. After
that, many Sparse DFTs with these input parameters can be computed (on
different input vectors). This concept is similar to FFTW's concept of plans.

You can create a plan with a call to ``sfft_plan``::

    sfft_plan* sfft_make_plan(int n, int k, sfft_version version,
                              int fftw_optimization);

The call returns a pointer to a struct of type ``sfft_plan``, which has to be
manually freed with ``sfft_free_plan``. Parameters of ``sfft_make_plan`` are:

``n``
    The size of the input vector.
``k``
    The number of frequencies in the signal, i.e. the signal's *sparsity*.
``version``
    The SFFT algorithm version to use. Either ``SFFT_VERSION_1``,
    ``SFFT_VERSION_2``, or ``SFFT_VERSION_3``.
``fftw_optimization``
    FFTW optimization level. Usually one of ``FFTW_MEASURE`` and 
    ``FFTW_ESTIMATE``. Since experiments showed that there is little benefit 
    in using the more expensive ``FFTW_MEASURE``, the best choice is typically
    ``FFTW_ESTIMATE``.


Creating Input Vectors
++++++++++++++++++++++

The storage for SFFT input vectors has to allocated using ``sfft_malloc``::

    void* sfft_malloc(size_t s);

The reason for this is that the implementation requires a specific memory
alignment on the input vectors. You can use ``sfft_malloc`` as a drop-in
replacement for ``malloc``.

Input vectors should be of type ``complex_t``, which is a typedef to the C
standard library's type ``double complex``.

Storage allocated with ``sfft_malloc`` must be freed with this function::

    void sfft_free(void*);


Creating the Output Datastructure
+++++++++++++++++++++++++++++++++

The output of the SFFT is stored in an associative array that maps frequency
coordinates to coefficients.  The array should be of type ``sfft_output``,
which is a typedef to an ``std::unordered_map``. Before executing the SFFT
plans, you need to create the output datastructure. A pointer to it is passed
to the SFFT execution call and the datastructure filled with the result.


Computing a Single Sparse DFT
+++++++++++++++++++++++++++++

Once a plan is created, input vectors are created filled with data, and an
output object was allocated, the SFFT plans can be executed. The function for
this is::

    void sfft_exec(sfft_plan* plan, complex_t* in, sfft_output* out);

Parameters should be self-explanatory. After execution of this function, the
output of the DFT is stored in ``*out``.


Computing Multiple Sparse DFTs
++++++++++++++++++++++++++++++

If you want to run multiple SFFT calls on different inputs (but with the same
input sizes), you can use ``sfft_exec_many`` to run the calls in parallel::

    void sfft_exec_many(sfft_plan* plan, 
                        int num, complex_t** in, sfft_output* out);

The function is very similar to ``sfft_exec``, but you can pass it put ``num``
input-vectors and ``num`` output-objects. The SFFT library used OpenMP for
parallelization; thus, you can use either the environment variable
``OMP_NUM_THREADS`` or OpenMP library functions to adjust the number of
threads. Be careful: do *not* use different thread number configuration for
the  call to ``sfft_make_plan`` and ``sfft_exec_many``. Otherwise your
program will crash!


SFFT Versions
-------------

Currently, three different SFFT versions are implemented: SFFT v1, v2, and v3.

SFFT v3 is the algorithm of choice when your input signals are exactly-sparse;
that is, there is no additional noise in the signals. SFFT v3 will not work
with noisy signals.


SFFT v1 and v2 can also be applied to noisy signals, but they only work with
certain input parameter combinations. Valid input parameters combinations:

============ ========
Signal Size  Sparsity
============ ========
8192         50
16384        50
32768        50
65536        50
131072       50
262144       50
524288       50
1048576      50
2097152      50
4194304      50
8388608      50
16777216     50
4194304      50
4194304      100
4194304      200
4194304      500
4194304      1000
4194304      2000
4194304      2500
4194304      4000
============ ========
