import pylab, os
from numpy import *
from numpy.linalg import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def read_coords(path):
    lines = file(path).readlines()

    cdim = int(lines[0].split()[1])
    sdim = int(lines[1].split()[1])

    lines = lines[2:]
    
    x = array([line.strip().split()[0] for line in lines], float).T
    y = array([line.strip().split()[1] for line in lines], float).T
    z = array([line.strip().split()[2] for line in lines], float).T
    nc = x.size

    if (nc != cdim*sdim):
        print 'Error: number of coordinates read'

    x = reshape(x, (cdim,sdim))
    y = reshape(y, (cdim,sdim))
    z = reshape(z, (cdim,sdim))
 
    return cdim, sdim, x, y, z

def read_mode(path):
    lines = file(path).readlines()

    cdim = int(lines[0].split()[1])
    sdim = int(lines[1].split()[1])

    lines = lines[2:]
    
    V = array([line.strip().split()[0] for line in lines], float).T
    nc = V.size

    if (nc != cdim*sdim):
        print 'Error: number of coordinates read'

    V = reshape(V, (cdim,sdim))
 
    return cdim, sdim, V

def split_blade(cdim, sdim, x, y, z):
    # assume that o-mesh is cut at LE or TE
    if mod(cdim,2) == 1:
        xu = x[:(cdim-1)/2+1,:]
        yu = y[:(cdim-1)/2+1,:]
        zu = z[:(cdim-1)/2+1,:]

        xl = x[(cdim-1)/2:,:]
        yl = y[(cdim-1)/2:,:]
        zl = z[(cdim-1)/2:,:]
    else:
        xu = x[:(cdim)/2+1,:]
        yu = y[:(cdim)/2+1,:]
        zu = z[:(cdim)/2+1,:]

        xl = x[(cdim)/2-1:,:]
        yl = y[(cdim)/2-1:,:]
        zl = z[(cdim)/2-1:,:]

    xl = xl[::-1,:]
    yl = yl[::-1,:]
    zl = zl[::-1,:]

    return xu, yu, zu, xl, yl, zl

def xyz2st(x,y,z):
    s0 = sum(sqrt((x[1:,0]-x[0:-1,0])**2 +\
                  (y[1:,0]-y[0:-1,0])**2 +\
                  (z[1:,0]-z[0:-1,0])**2))
    sN = sum(sqrt((x[1:,-1]-x[0:-1,-1])**2 +\
                  (y[1:,-1]-y[0:-1,-1])**2 +\
                  (z[1:,-1]-z[0:-1,-1])**2))

    if (s0 > sN):
        s = cumsum(sqrt((x[1:,0]-x[0:-1,0])**2 +\
                        (y[1:,0]-y[0:-1,0])**2 +\
                        (z[1:,0]-z[0:-1,0])**2))
        s = hstack((0.,s))
    else: 
        s = cumsum(sqrt((x[1:,-1]-x[0:-1,-1])**2 +\
                        (y[1:,-1]-y[0:-1,-1])**2 +\
                        (z[1:,-1]-z[0:-1,-1])**2))
        s = hstack((0.,s))

    t0 = sum(sqrt((x[0,1:]-x[0,0:-1])**2 +\
                  (y[0,1:]-y[0,0:-1])**2 +\
                  (z[0,1:]-z[0,0:-1])**2))
    tN = sum(sqrt((x[-1,1:]-x[-1,0:-1])**2 +\
                  (y[-1,1:]-y[-1,0:-1])**2 +\
                  (z[-1,1:]-z[-1,0:-1])**2))

    if (t0 > tN):
        t = cumsum(sqrt((x[0,1:]-x[0,0:-1])**2 +\
                        (y[0,1:]-y[0,0:-1])**2 +\
                        (z[0,1:]-z[0,0:-1])**2))
        t = hstack((0.,t))
    else: 
        t = cumsum(sqrt((x[-1,1:]-x[-1,0:-1])**2 +\
                        (y[-1,1:]-y[-1,0:-1])**2 +\
                        (z[-1,1:]-z[-1,0:-1])**2))
        t = hstack((0.,t))

    return s, t

def calcNormals(x, y, z):
    nr = x.shape[0]-1
    nc = x.shape[1]-1
    # face areas
    areas = zeros((nr,nc))
    v1 = zeros((nr,nc,3))
    v1[:,:,0] = x[1:,:-1] - x[:-1,:-1]
    v1[:,:,1] = y[1:,:-1] - y[:-1,:-1]
    v1[:,:,2] = z[1:,:-1] - z[:-1,:-1]
    v2 = zeros((nr,nc,3))
    v2[:,:,0] = x[:-1,1:] - x[:-1,:-1]
    v2[:,:,1] = y[:-1,1:] - y[:-1,:-1]
    v2[:,:,2] = z[:-1,1:] - z[:-1,:-1]
    n = 0.5*cross(v1.reshape(nr*nc,3),v2.reshape(nr*nc,3)).reshape(nr,nc,3)
    areas += sqrt(n[:,:,0]**2 + n[:,:,1]**2 + n[:,:,2]**2)
    v1[:,:,0] = x[1:,1:] - x[:-1,1:]
    v1[:,:,1] = y[1:,1:] - y[:-1,1:]
    v1[:,:,2] = z[1:,1:] - z[:-1,1:]
    v2[:,:,0] = x[1:,1:] - x[1:,:-1]
    v2[:,:,1] = y[1:,1:] - y[1:,:-1]
    v2[:,:,2] = z[1:,1:] - z[1:,:-1]
    n = 0.5*cross(v1.reshape(nr*nc,3),v2.reshape(nr*nc,3)).reshape(nr,nc,3)
    areas += sqrt(n[:,:,0]**2 + n[:,:,1]**2 + n[:,:,2]**2)

    areas = column_stack((areas[:,0],areas,areas[:,-1]))
    areas = row_stack((areas[0,:],areas,areas[-1,:]))

    # face normals - construct with face diagonals
    v1[:,:,0] = x[1:,1:] - x[:-1,:-1]
    v1[:,:,1] = y[1:,1:] - y[:-1,:-1]
    v1[:,:,2] = z[1:,1:] - z[:-1,:-1]
    v2[:,:,0] = x[1:,:-1] - x[:-1,1:]
    v2[:,:,1] = y[1:,:-1] - y[:-1,1:]
    v2[:,:,2] = z[1:,:-1] - z[:-1,1:]
    nf = 0.5*cross(v1.reshape(nr*nc,3),v2.reshape(nr*nc,3)).reshape(nr,nc,3)
    nf = nf/sqrt(nf[:,:,0]**2 + nf[:,:,1]**2 + nf[:,:,2]**2).reshape(nr,nc,1).repeat(3,2)

    nfx = column_stack((nf[:,0,0],nf[:,:,0],nf[:,-1,0]))
    nfy = column_stack((nf[:,0,1],nf[:,:,1],nf[:,-1,1]))
    nfz = column_stack((nf[:,0,2],nf[:,:,2],nf[:,-1,2]))
    nfx = row_stack((nfx[0,:],nfx,nfx[-1,:]))
    nfy = row_stack((nfy[0,:],nfy,nfy[-1,:]))
    nfz = row_stack((nfz[0,:],nfz,nfz[-1,:]))

    n = zeros((nr+1,nc+1,3))
    n[:,:,0] = 0.25*(nfx[:-1,:-1]*areas[:-1,:-1] + nfx[1:,1:]*areas[1:,1:] + 
                     nfx[:-1,1:]*areas[:-1,1:] + nfx[1:,:-1]*areas[1:,:-1])
    n[:,:,1] = 0.25*(nfy[:-1,:-1]*areas[:-1,:-1] + nfy[1:,1:]*areas[1:,1:] + 
                     nfy[:-1,1:]*areas[:-1,1:] + nfy[1:,:-1]*areas[1:,:-1])
    n[:,:,2] = 0.25*(nfz[:-1,:-1]*areas[:-1,:-1] + nfz[1:,1:]*areas[1:,1:] + 
                     nfz[:-1,1:]*areas[:-1,1:] + nfz[1:,:-1]*areas[1:,:-1])

    n = n/sqrt(n[:,:,0]**2 + n[:,:,1]**2 + n[:,:,2]**2).reshape(nr+1,nc+1,1).repeat(3,2)

    # make sure periodicity is correct
    n[-1,:,:] = n[0,:,:]

    return n

def calcNormalsCamber(x, y, z):
    nr = x.shape[0]-1
    nc = x.shape[1]-1
    # face areas
    areas = zeros((nr,nc))
    v1 = zeros((nr,nc,3))
    v1[:,:,0] = x[1:,:-1] - x[:-1,:-1]
    v1[:,:,1] = y[1:,:-1] - y[:-1,:-1]
    v1[:,:,2] = z[1:,:-1] - z[:-1,:-1]
    v2 = zeros((nr,nc,3))
    v2[:,:,0] = x[:-1,1:] - x[:-1,:-1]
    v2[:,:,1] = y[:-1,1:] - y[:-1,:-1]
    v2[:,:,2] = z[:-1,1:] - z[:-1,:-1]
    n = 0.5*cross(v1.reshape(nr*nc,3),v2.reshape(nr*nc,3)).reshape(nr,nc,3)
    areas += sqrt(n[:,:,0]**2 + n[:,:,1]**2 + n[:,:,2]**2)
    v1[:,:,0] = x[1:,1:] - x[:-1,1:]
    v1[:,:,1] = y[1:,1:] - y[:-1,1:]
    v1[:,:,2] = z[1:,1:] - z[:-1,1:]
    v2[:,:,0] = x[1:,1:] - x[1:,:-1]
    v2[:,:,1] = y[1:,1:] - y[1:,:-1]
    v2[:,:,2] = z[1:,1:] - z[1:,:-1]
    n = 0.5*cross(v1.reshape(nr*nc,3),v2.reshape(nr*nc,3)).reshape(nr,nc,3)
    areas += sqrt(n[:,:,0]**2 + n[:,:,1]**2 + n[:,:,2]**2)

    areas = column_stack((areas[:,0],areas,areas[:,-1]))
    areas = row_stack((areas[0,:],areas,areas[-1,:]))

    # face normals - construct with face diagonals
    v1[:,:,0] = x[1:,1:] - x[:-1,:-1]
    v1[:,:,1] = y[1:,1:] - y[:-1,:-1]
    v1[:,:,2] = z[1:,1:] - z[:-1,:-1]
    v2[:,:,0] = x[1:,:-1] - x[:-1,1:]
    v2[:,:,1] = y[1:,:-1] - y[:-1,1:]
    v2[:,:,2] = z[1:,:-1] - z[:-1,1:]
    nf = 0.5*cross(v1.reshape(nr*nc,3),v2.reshape(nr*nc,3)).reshape(nr,nc,3)
    nf = nf/sqrt(nf[:,:,0]**2 + nf[:,:,1]**2 + nf[:,:,2]**2).reshape(nr,nc,1).repeat(3,2)

    nfx = column_stack((nf[:,0,0],nf[:,:,0],nf[:,-1,0]))
    nfy = column_stack((nf[:,0,1],nf[:,:,1],nf[:,-1,1]))
    nfz = column_stack((nf[:,0,2],nf[:,:,2],nf[:,-1,2]))
    nfx = row_stack((nfx[0,:],nfx,nfx[-1,:]))
    nfy = row_stack((nfy[0,:],nfy,nfy[-1,:]))
    nfz = row_stack((nfz[0,:],nfz,nfz[-1,:]))

    n = zeros((nr+1,nc+1,3))
    n[:,:,0] = 0.25*(nfx[:-1,:-1]*areas[:-1,:-1] + nfx[1:,1:]*areas[1:,1:] + 
                     nfx[:-1,1:]*areas[:-1,1:] + nfx[1:,:-1]*areas[1:,:-1])
    n[:,:,1] = 0.25*(nfy[:-1,:-1]*areas[:-1,:-1] + nfy[1:,1:]*areas[1:,1:] + 
                     nfy[:-1,1:]*areas[:-1,1:] + nfy[1:,:-1]*areas[1:,:-1])
    n[:,:,2] = 0.25*(nfz[:-1,:-1]*areas[:-1,:-1] + nfz[1:,1:]*areas[1:,1:] + 
                     nfz[:-1,1:]*areas[:-1,1:] + nfz[1:,:-1]*areas[1:,:-1])

    n = n/sqrt(n[:,:,0]**2 + n[:,:,1]**2 + n[:,:,2]**2).reshape(nr+1,nc+1,1).repeat(3,2)

    return n
