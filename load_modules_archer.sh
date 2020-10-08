#!/bin/bash

module swap PrgEnv-cray PrgEnv-gnu
module load cmake/3.16.0
module load python-compute/3.6.0_gcc6.1.0
module swap gcc/6.1.0 gcc/7.3.0
module load boost/1.60
