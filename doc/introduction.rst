Introduction
============

The *Sparse Fast Fourier Transform* is a DFT algorithm specifically designed
for signals with a sparse frequency domain. This library is a high-performance
C++ implementation of versions 1, 2, and 3 of the different SFFT variants. 


When Should I use the SFFT library?
-----------------------------------

You should use the SFFT library when you want to compute the `Discrete Fourier
Transform`_ of a signal and only a few frequency components
occur in the signal. Your signal may be noisy or not, but currently there are
some limitations for noisy signals (see :ref:`current-state`).

.. _Discrete Fourier Transform: 
      http://en.wikipedia.org/wiki/Discrete_Fourier_transform



Target Platform
---------------

The SFFT library was optimized to run on modern x86 desktop CPUs with SSE
support. Optionally the implementation can use the Intel IPP library, which is
only available on Intel platforms.

.. _current-state:

Limitations and Known Bugs
--------------------------

The SFFT library features implementations of SFFT v1, v2, and v3. SFFT v1 and
v2 currently only work with a few specific input parameters. SFFT v3 cannot
handle signals with noise. 

There are no known bugs so far.

Disclaimer
----------

The current SFFT implementation is in an experimental state. It is NOT
intended to be used as a drop-in replacement for the FFT library of your choice.
Be prepared to find bugs. There is absolutely NO WARRANTY for the correct
functioning of this software.


Credits
-------

The original SFFT sourcecode was developed by Haitham Hassanieh, Piotr Indyk,
Dina Katabi, and Eric Price at the Computer Science and Artifical Intelligence
Lab at MIT. The original sourcecode and contact information can be found at
their website `Sparse Fast Fourier Transform Website`_.

Performance optimizations were developed by Jörn Schumacher
as part of his `Master Thesis Project`_ at the
Computer Science Department of ETH Zurich in 2013, under the supervision of
Prof. `Markus Püschel`_.

Contact Information
-------------------

If you are interested in the theory behind the Sparse Fast Fourier Transform,
contact the inventors of the SFFT, Haitham Hassanieh, Piotr Indyk, Dina Katabi,
and Eric Price, at their `Sparse Fast Fourier Transform Website`_.

If you are interested in performance optimizations that were applied, contact
Jörn Schumacher at joerns@student.ethz.ch.

.. _Sparse Fast Fourier Transform Website:
        http://groups.csail.mit.edu/netmit/sFFT/
.. _Master Thesis Project: http://www.spiral.net/software/sfft.html
.. _Markus Püschel: http://www.inf.ethz.ch/personal/markusp/

