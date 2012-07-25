from numpy import *

# write the blade surface in tecplot format
def write_blade_surf(x,y,z,fname):

    cdim = x.shape[0]
    sdim = x.shape[1]
    np = cdim*sdim
    ne = 2*(cdim-1)*(sdim-1)

    f = open(fname,'w')
    f.write('VARIABLES = "X", "Y", "Z", ')
    f.write('\n')
    f.write('ZONE N='+str(np)+', E='+str(ne)+', DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n')
    
    # write out the data
    for i in range(cdim):
        for j in range(sdim):
            f.write('%20.8E'*3 % (x[i,j],y[i,j],z[i,j]))
            f.write('\n')
    
    # write out the elements
    for i in range(cdim-1):
        for j in range(sdim-1):
            a = sdim*i + (j+1)
            b = sdim*i + (j+2)
            c = sdim*(i+1) + (j+1)
            d = sdim*(i+1) + (j+2)
            f.write('%d  '*3 % (a,b,c))
            f.write('\n')
            f.write('%d  '*3 % (c,b,d))
            f.write('\n')

    f.close()
