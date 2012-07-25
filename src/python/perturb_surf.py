import read_blade, write_tecplot
from numpy import *
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.interpolate import *

def perturb_chev(rpath, wpath, M):

    # M - amplitudes of the Chevyshev modes

    cdim, sdim, x, y, z = read_blade.read_coords(rpath)

    # split the blade
    xu,yu,zu,xl,yl,zl = read_blade.split_blade(cdim,sdim,x,y,z)

    # compute the s coordinates
    tck, su = splprep([xu[:,0],yu[:,0]],s=0)
    tck, sl = splprep([xl[:,0],yl[:,0]],s=0)

    # form the modes
    nmodes = len(M)
    cheb = zeros((nmodes,x.shape[0],x.shape[1]))
    for i in arange(nmodes):
        n = i + 1
        if mod(n,2) == 0:
            gu =  (1-2*su-cos((n+1)*arccos(1-2*su)))/(n+1.)
            gl = -(1-2*sl-cos((n+1)*arccos(1-2*sl)))/(n+1.)
            cheb[i,:,:] = tile(hstack((gu,gl[::-1][1:])),(sdim,1)).T
        else:
            gu =  (1-cos((n+1)*arccos(1-2*su)))/(n+1.)
            gl = -(1-cos((n+1)*arccos(1-2*sl)))/(n+1.)
            cheb[i,:,:] = tile(hstack((gu,gl[::-1][1:])),(sdim,1)).T

    np = read_blade.calcNormals(x,y,z)

    for i in arange(nmodes):
        xp = x + np[:,:,0]*cheb[i,:,:]*M[i]
        yp = y + np[:,:,1]*cheb[i,:,:]*M[i]
        zp = z + np[:,:,2]*cheb[i,:,:]*M[i]

    pylab.plot(x[:,0],y[:,0])
    pylab.plot(xp[:,0],yp[:,0])
    pylab.axis('Equal')
    pylab.show()
     
    # write out the blade surface
    f = open(wpath,'w')
    f.write('CDIM:       %d\n' % cdim)
    f.write('SDIM:       %d\n' % sdim)
    for i in arange(cdim):
        for j in arange(sdim):
            f.write('%20.8E' * 3 % (xp[i,j],yp[i,j],zp[i,j]))
            f.write('\n')
    
    f.close()

    # write out to tecplot format
    write_tecplot.write_blade_surf(xp,yp,zp,'tec_blade.dat')

def perturb_rot(rpath, wpath, dth):

    # dht - amount to rotate blade

    cdim, sdim, x, y, z = read_blade.read_coords(rpath)

    # split the blade
    xu,yu,zu,xl,yl,zl = read_blade.split_blade(cdim,sdim,x,y,z)

    xle = xu[0,0]; yle = yu[0,0]
    xte = xu[-1,0]; yte = yu[-1,0]
    th = arctan2(y-yle,x-xle)
    r = sqrt((x-xle)**2 + (y-yle)**2)
    xp = xle + r*cos(th+dth)
    yp = yle + r*sin(th+dth)
    zp = copy(z)

    pylab.plot(x[:,0],y[:,0])
    pylab.plot(xp[:,0],yp[:,0])
    pylab.axis('Equal')
    pylab.show()
     
    # write out the blade surface
    f = open(wpath,'w')
    f.write('CDIM:       %d\n' % cdim)
    f.write('SDIM:       %d\n' % sdim)
    for i in arange(cdim):
        for j in arange(sdim):
            f.write('%20.8E' * 3 % (xp[i,j],yp[i,j],zp[i,j]))
            f.write('\n')
    
    f.close()

    # write out to tecplot format
    write_tecplot.write_blade_surf(xp,yp,zp,'tec_blade.dat')
