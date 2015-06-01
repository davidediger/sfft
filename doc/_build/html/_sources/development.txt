.. _development:

Development
===========

Development and Benchmark Tools
-------------------------------

The SFFT library includes some useful tools for development and benchmarking.
To enable them, you have to configure with the ``--develop`` flag. Then, the
following programs will be built additionally:


``sfft-cachemisses``
    Runs an SFFT on random input. The tool is handy when used with Valgrind's
    cachegrind tool. The program includes some instructions to disable 
    valgrind during the input-generation and planning phases. Thus, when the
    program is analyzed with cachegrind, only the execution phase will be
    observed.
``sfft-instruction_count``
	Counts the floating point instructions of the specified SFFT call
	(configured with program parameters, see below) and prints them. When the
	configuration option ``--profile`` was defined, this will also print a
	profile	of the SFFT call.
``sfft-profiling``
	Another program that runs a configurable SFFT call. This program will be
	compiled with the profiling flags ``-pg``, so that it can be analyzed with
	the ``gprof`` profiling tool.
``sfft-timing``
	A program that accurately measures the runtime of the specified SFFT call.
	This can be used by benchmark scripts.
``sfft-timing_many``
	Similar to ``sfft-timing``, but measures the parallel execution of
	multiple SFFT calls.
``sfft-verification``
	This program runs the specified SFFT call and checks that the output is
	correct. This is useful for testing.


All of the programs run one or many SFFT executions. Random input data is
generated automatically. The programs share the following common options:

``-n SIZE``
	The size of the input signal.
``-k NUMBER``
	Number of frequencies generated in the random input signal.
``-r REPETITIONS``
	*NOT available for sfft-timing_many.* Allows to compute multiple SFFTs.
	Default: 1. .
``-i NUM``
	*Only available for sfft-timing_many.* Generate NUM inputs. 
``-s``
	*Only available for sfft-timing_many.* Do not share data between
	threads. This is slower.
``-v VERSION``
	Selects the algorithm version to use. ``VERSION`` is either 1, 2, or 3. ``
``-o``
	When ``-o`` is used, ``FFTW_MEASURE`` is used for FFTW calls instead of
	``FFTW_ESTIMATE``.
``-h``
	Displays help.


An Overview of the Sourcecode
-----------------------------

Here is an overview of the purpose of different sourcefiles:

*cachemisses.cc, timing.cc, timing_many.cc, instruction_count.cc, verification.cc, simulation.[cc,h]*
	The ``main`` routines and some support code for all development tools
	are located in these files.
*computefourier-1.0-2.0.[cc,h]*
	Algorithm sourcecode for SFFT v1 and v2.
*computefourier-3.0.[cc,h]*
	Algorithm sourcecode for SFFT v3.
*fft.h, common.[cc,h], utils.[cc,h]*
	Some common code and datatypes.
*fftw.[cc,h]*
	Interface code for FFTW calls.
*filters.[cc,h]*
	The routines to generate filter vectors are in here.
*intrinsics.h*
	Some compiler-specific abstractions to include the correct intrinsics
	header.
*parameters.[cc,h]*
	Parameter configuration for SFFT v1, v2.
*profiling_tools.h*
	Some preprocessor tools to allow profiling, used when compiled with
	``--profile``.
*roofline.cc*
	A program to use with the roofline tool *perfplot*. Can be built with
	``tools/build-roofline.sh``.
*sfft.[cc,h]*
	User interface code and basic datastructures. The headerfile is to be
	included by users.
*timer.[cc,h]*
	Functions for accurate timing, used by ``sfft-timing``.
*flopcount/*
	Files in this directory are used to count floating point operations, used
	by ``sfft-instruction_count``.

