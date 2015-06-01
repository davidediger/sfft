#! /usr/bin/env python
# encoding: utf-8

import os

# the following two variables are used by the target "waf dist"
VERSION='0.1.0'
APPNAME='sfft'

OPTIMIZATION_FLAGS = ['-O3', '-ffast-math', '-march=native', '-fopenmp']
OPENMP_FLAGS = '-fopenmp'
WARNING_FLAGS = ['-Wall']
DEBUG_FLAGS = ['-g']
PROFILING_FLAGS = ['-pg']

# these variables are mandatory ('/' are converted automatically)
top = '.'
out = 'build'

def options(opt):
  opt.load('compiler_cxx')
	
  opt.add_option('--develop', action='store_true', default=False, help='build development and benchmark tools')
  opt.add_option('--debug', action='store_true', default=False, help='compile in debug mode')
  opt.add_option('--profile', action='store_true', default=False, help='add source-level profiling to instruction counting programs')
  opt.add_option('--without-ipp', action='store_true', default=False, help='do not the Intel Performance Primitives library')

def configure(conf):
  conf.load('compiler_cxx')
  conf.check_cfg(package='valgrind', args='--cflags')
  conf.check_cfg(package='fftw3', args=['--cflags', '--libs'], uselib_store="FFTW3")
  conf.check_cxx(lib='m', uselib_store='M') 
  conf.check_cxx(lib='rt', uselib_store='RT')
  conf.check_cxx(lib='gomp', uselib_store='OPENMP')
  if not conf.options.without_ipp:
    conf.check_cxx(lib= ['ippvm', 'ipps', 'pthread'], uselib_store="IPPVM", define_name="HAVE_IPP")

  conf.check_cxx(cxxflags=OPTIMIZATION_FLAGS, cflags=OPTIMIZATION_FLAGS,
                 linkflags=OPTIMIZATION_FLAGS, uselib_store="OPTIMIZATION")
  conf.check_cxx(cxxflags=WARNING_FLAGS, cflags=WARNING_FLAGS,
                 linkflags=WARNING_FLAGS, uselib_store="WARNINGS")
  conf.check_cxx(cxxflags=PROFILING_FLAGS, cflags=PROFILING_FLAGS,
                 linkflags=PROFILING_FLAGS, uselib_store="PROFILING")

  if conf.options.develop:
    conf.env.develop = True

  if conf.options.debug:
    conf.check_cxx(cxxflags=DEBUG_FLAGS, cflags=DEBUG_FLAGS,
                   linkflags=DEBUG_FLAGS, uselib_store="DEBUG")

  if conf.options.profile:
    conf.env.profile = True

	
def build(bld):
  exclude_files =  ['src/cachemisses.cc', 'src/instruction_count.cc', 'src/timing.cc', 'src/timing_many.cc', 'src/verification.cc', 'src/roofline.cc']
  includes = ['src/flopcount/']

  common_use = "IPPVM OPENMP FFTW3 FFTW3F M RT OPTIMIZATION WARNINGS DEBUG"
	
  source_files = bld.path.ant_glob('src/*.cc', excl=exclude_files)
  bld.objects(source=source_files, includes=includes, target='sfft-objects', use=common_use)

  bld.stlib(source=source_files, target='sfft', includes=includes, features='cxx', use=common_use)
  bld.shlib(source=source_files, target='sfft', includes=includes, features='cxx', use=common_use)

  if bld.env.develop:
    for program in ['timing', 'timing_many', 'verification', 'cachemisses']:
      bld.program(features='cxx cxxprogram',
                  source='src/'+program+'.cc',
                  target='sfft-'+program,
                  includes=includes,
                  use='sfft-objects ' + common_use)

    flopcount_defines = ["USE_FLOPCOUNTER"]
    if bld.env.profile:
      flopcount_defines += ["USE_PROFILING"]
    bld.objects(source=source_files+['src/flopcount/flopcount.cpp'], includes=["src/flopcount/"], 
                defines=flopcount_defines, target='sfft-flopcount-objects', use=common_use)
    bld.program(features='cxx cxxprogram',
                source='src/instruction_count.cc',
                target='sfft-instruction_count',
                includes=includes,
                defines='USE_FLOPCOUNTER',
                use='sfft-flopcount-objects ' + common_use)

    bld.objects(source=source_files, includes=includes,
                target='sfft-profile-objects', use='PROFILING '+common_use)
    bld.program(features='cxx cxxprogram',
                source='src/timing.cc',
                target='sfft-profiling',
                includes=includes,
                use='sfft-profile-objects ' + 'PROFILING ' + common_use)
    
