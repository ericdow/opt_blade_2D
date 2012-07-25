import perturb_surf
import os, shutil
from filelock import FileLock
from numpy import *

def modify(src,inp,rundir,cg_mesh_orig,cg_mesh_mod,rpath,wpath,M):
    # src          : location of perturb_mesh fortran code
    # inp          : location of CGNS files
    # rundir       : current running directory
    # cg_mesh_orig : name of original CGNS file
    # cg_mesh_mod  : name of modified CGNS file
    # rpath        : name of original blade surface file
    # wpath        : name of modified blade surface file
    # M            : perturbation mode amplitudes

    with FileLock(inp+cg_mesh_orig):
        # copy CGNS to rundir
        shutil.copy(inp+cg_mesh_orig,rundir)

    os.chdir(rundir)

    # convert CGNS files to HDF format
    os.system('adf2hdf ' + cg_mesh_orig)
    
    # copy CGNS file
    shutil.copy(cg_mesh_orig, cg_mesh_mod)
    
    # write mesh surface out
    os.system(src+'perturb_mesh_2D '+cg_mesh_orig)
    
    # perturb the mesh surface
    # perturb_surf.perturb_chev(rundir+rpath, rundir+wpath, M)
    perturb_surf.perturb_rot(rundir+rpath, rundir+wpath, 1*pi/180.)
    
    # read in perturbation to CGNS mesh
    os.system(src+'perturb_mesh_2D '+cg_mesh_orig+' '+cg_mesh_mod)
    
    # convert CGNS file back to ADF format
    os.system('hdf2adf ' + cg_mesh_mod)

    # remove the original CGNS file
    os.remove(cg_mesh_orig)
