import mod_mesh
import os, shutil
from numpy import *

# ut_src = '/mnt/pwfiles/ericdow/utcfd/bin/'
src = '/home/ericdow/code/opt_blade_2D/src/'
inp = '/home/ericdow/code/opt_blade_2D/input/'
rundir = '/home/ericdow/code/opt_blade_2D/input/tmp/'

cg_mesh_orig = 'sc10.cgns'
cg_mesh_mod  = 'sc10_mod.cgns'

rpath = 'blade_surf.dat'
wpath = 'blade_surf_mod.dat'

# number of iterations to run each mesh
niter = 3000
# number of cores to run on
np = 1

# generate modified mesh
nmodes = 5
M = [.005,.005,.005,.005,.005]
mod_mesh.modify(src,inp,rundir,cg_mesh_orig,cg_mesh_mod,rpath,wpath,M)

