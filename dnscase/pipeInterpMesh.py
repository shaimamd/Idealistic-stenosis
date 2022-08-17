import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math as mt
import struct
import sys


def pointsinfda():
    ## Simulation parameters
    D = 1
    Rmax = D / 2
    # Radial points
    nR = 10  # Number of Radii
    R = np.linspace(0, Rmax, nR)

    # Azimuthal points
    nth = 15
    arclenght = (R * 3 * mt.pi) - (R * 2 * mt.pi) / nth
    # lastv = arclenght[-1] / Rmax
    # th = np.linspace(0, lastv, nth)
    th = np.linspace(0, 2 * mt.pi, nth)

    # Streamwise points for long expansion

    zn = np.arange(0,54,1)
    nz = len(zn)

    print("We kept", nz, "sections in Z")

    ##############################################################

    nums = len(th) * len(R) * len(zn)

    print("We kept totally", nums, "points in the domain")

    xx = np.zeros((nR, nth))
    yy = np.zeros((nR, nth))
    theta = np.zeros((nR, nth))

    for i in range(nR):
        for j in range(nth):
            xx[i, j] = R[i] * mt.cos(th[j])
            yy[i, j] = R[i] * mt.sin(th[j])
            theta[i, j] = th[j]

    # flatten in polar arrangement
    x_pts = xx.flatten('F')  # xx(:)
    y_pts = yy.flatten('F')  # yy(:)
    th_pts = theta.flatten('F')  # theta(:)

    # Copy paste same xy points to every z-section
    repetitions = nz
    xx2 = np.tile(x_pts, (repetitions, 1))
    yy2 = np.tile(y_pts, (repetitions, 1))
    angle2 = np.tile(th_pts, (repetitions, 1))
    z = np.repeat(zn, int(nR * nth))

    x2 = xx2.flatten('C')
    y2 = yy2.flatten('C')
    angle = angle2.flatten('C')

    # check lengths:
    print("length of x, y and z is", len(x2), len(y2), len(z))

    # Modify shape for the nozzle radius
    # nozzle starts at -15.67 in the original
    # throat starts at -10 in the original
    zlong1 = 2.0
    zlong2 =4.0

    f0 = 1 / 4
    L=2.0
    Xo=2.0
    x = np.zeros(nums)
    y = np.zeros(nums)

    for i in range(len(z)):
        if z[i] <= zlong2 and z[i] >= zlong1:  # nozzle section
            s=(1-f0*(1-mt.cos(2*mt.pi*((z[i])-Xo)/(L))))
            e=f0*(1+mt.cos(2*mt.pi*((z[i])-Xo)/(L)))/10
           
            x[i] = x2[i] *s
            y[i] = y2[i] *s+e

        else:  # remaining places
            x[i] = x2[i]
            y[i] = y2[i]

    return x, y, z, angle, nums


# ---------------------------------------------------------------------------------
# Now to write these points into file :

class point:
    """class defining point variables"""

    def __init__(self, ldim):
        self.pos = np.zeros((ldim))


class pset:
    """class containing data of the point collection"""

    def __init__(self, ldim, npoints):
        self.ldim = ldim
        self.npoints = npoints
        self.pset = [point(ldim) for il in range(npoints)]


def set_pnt_pos(data, il, lpos):
    """set position of the single point"""
    lptn = data.pset[il]
    data_pos = getattr(lptn, 'pos')
    for jl in range(data.ldim):
        data_pos[jl] = lpos[jl]


def write_int_pos(fname, wdsize, emode, data):
    """ write point positions to the file"""
    # open file
    outfile = open(fname, 'wb')

    # word size
    if (wdsize == 4):
        realtype = 'f'
    elif (wdsize == 8):
        realtype = 'd'

    # header
    header = '#iv1 %1i %1i %10i ' % (wdsize, data.ldim, data.npoints)
    header = header.ljust(32)
    outfile.write(header.encode('utf-8'))

    # write tag (to specify endianness)
    # etagb = struct.pack(emode+'f', 6.54321)
    # outfile.write(etagb)
    outfile.write(struct.pack(emode + 'f', 6.54321))

    # write point positions
    for il in range(data.npoints):
        lptn = data.pset[il]
        data_pos = getattr(lptn, 'pos')
        outfile.write(struct.pack(emode + data.ldim * realtype, *data_pos))

    # ----------------------------------------------------------------------


if __name__ == "__main__":
    # initialise variables
    fname = 'int_pos'
    wdsize = 8
    # little endian
    emode = '<'
    # big endian
    # emode = '<'

    #
    # generate polar mesh over a circular disk
    xx, yy, zz, thetas, npts = pointsinfda()

    # allocate space
    ldim = 3  # 3d
    npoints = npts
    data = pset(ldim, npoints)
    print('Allocated {0} points'.format(npoints))

    # initialise point positions
    lpos = np.zeros(data.ldim)

    # set the positions in data
    il = 0

    for il in range(0, npts):
        lpos[0] = xx[il]
        lpos[1] = yy[il]
        lpos[2] = zz[il]
        set_pnt_pos(data, il, lpos)

    # write points to the file
    write_int_pos(fname, wdsize, emode, data)

    # plot the mesh
    plt.plot(xx, yy, '.b')
    plt.show()

    print("I think we're done")

# -----------------------------------------------------------------------

