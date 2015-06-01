#!/bin/bash

# This script calls cachegrind (a valgrind tool) on the given program,
# but the collection is only performed within the function
#      
#             Simulation::run()
#
# This way, it is suitable for profiling sFFT programs


valgrind --tool=cachegrind $*

