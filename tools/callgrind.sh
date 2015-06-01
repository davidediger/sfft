#!/bin/bash

# This script calls callgrind (a valgrind tool) on the given program,
# but the collection is only performed within the function
#      
#             Simulation::run()
#
# This way, it is suitable for profiling sFFT programs


valgrind --tool=callgrind --collect-atstart=no --toggle-collect="simulation::run()" $*

