 # Example script for reading int_fld data
# Kept similar to pymech
import struct
import numpy as np

class point:
    """class defining point variables"""

    def __init__(self,ldim,nstat):
        self.pos = np.zeros((ldim))
        self.stat = np.zeros((nstat))

class pset:
    """class containing data of the point collection"""

    def __init__(self,ldim,nstat,npoints):
        self.ldim = ldim
        self.nstat = nstat
        self.npoints = npoints
        self.re = []
        self.bsize = np.zeros((3))
        self.belem = []
        self.pord = []
        self.start_time = []
        self.end_time = []
        self.effav_time = []
        self.dt = []
        self.nrec = []
        self.int_time = []
        self.pset = [point(ldim,nstat) for il in range(npoints)]

def set_pnt_pos(data,il,lpos):
    """set position of the single point"""
    lptn = data.pset[il]
    data_pos = getattr(lptn,'pos')
    for jl in range(data.ldim):
            data_pos[jl] =  lpos[jl]

def read_int(infile,emode,nvar):
    """read integer array"""
    isize = 4
    llist = infile.read(isize*nvar)
    llist = list(struct.unpack(emode+nvar*'i', llist))
    return llist
def write_int_pos(fname,wdsize,emode,data):
    """ write point positions to the file"""
    # open file
    outfile = open(fname, 'wb')

    # word size
    if (wdsize == 4):
        realtype = 'f'
    elif (wdsize == 8):
        realtype = 'd'

    # header
    header = '#iv1 %1i %1i %10i ' %(wdsize,data.ldim,data.npoints)
    header = header.ljust(32)
    outfile.write(header.encode('utf-8'))

    # write tag (to specify endianness)
    #etagb = struct.pack(emode+'f', 6.54321)
    #outfile.write(etagb)
    outfile.write(struct.pack(emode+'f', 6.54321))

    #write point positions
    for il in range(data.npoints):
        lptn = data.pset[il]
        data_pos = getattr(lptn,'pos')
        outfile.write(struct.pack(emode+data.ldim*realtype, *data_pos))
def read_flt(infile,emode,wdsize,nvar):
    """read real array"""
    if (wdsize == 4):
        realtype = 'f'
    elif (wdsize == 8):
        realtype = 'd'
    llist = infile.read(wdsize*nvar)
    llist = np.frombuffer(llist, dtype=emode+realtype, count=nvar)
    return llist

def read_int_fld(fname):
    """read data from interpolation file"""
    # open file
    infile = open(fname, 'rb')
    # read header
    header = infile.read(460).split()

    # extract word size
    wdsize = int(header[1])

    # identify endian encoding
    etagb = infile.read(4)
    etagL = struct.unpack('<f', etagb)[0]; etagL = int(etagL*1e5)/1e5
    etagB = struct.unpack('>f', etagb)[0]; etagB = int(etagB*1e5)/1e5
    if (etagL == 6.54321):
        emode = '<'
    elif (etagB == 6.54321):
        emode = '>'

    # get simulation parameters
    ren = read_flt(infile,emode,wdsize,1)
    bsize = read_flt(infile,emode,wdsize,3)
    belem = read_int(infile,emode,3)
    pord = read_int(infile,emode,3)
    nstat = read_int(infile,emode,1)
    stime = read_flt(infile,emode,wdsize,4)
    nrec = read_int(infile,emode,1)
    itime = read_flt(infile,emode,wdsize,1)
    npoints = read_int(infile,emode,1)

    # create main data structure
    ldim = 3
    data = pset(ldim,nstat[0],npoints[0])

    # fill simulation parameters
    data.re = ren[0]
    data.bsize = bsize
    data.belem = belem
    data.pord = pord
    data.start_time = stime[0]
    data.end_time = stime[1]
    data.effav_time = stime[2]
    data.dt = stime[3]
    data.nrec = nrec[0]
    data.int_time = itime[0]

    # read coordinates
    for il in range(data.npoints):
        lpos = read_flt(infile,emode,wdsize,data.ldim)
        lptn = data.pset[il]
        data_pos = getattr(lptn,'pos')
        for jl in range(data.ldim):
            data_pos[jl] =  lpos[jl]

    # read statistics
    for il in range(data.nstat):
        for jl in range(data.npoints):
            lsts = read_flt(infile,emode,wdsize,1)
            lptn = data.pset[jl]
            data_sts = getattr(lptn,'stat')
            data_sts[il] =  lsts[0]

    return data
if __name__ == "__main__":
    # initialise variables
    fname = 'int_pos'
    wdsize = 8
    # little endian
    emode = '<'
    # big endian
    #emode = '<'

   
    
    nz = np.arange(0,54,1)
    # Radial points
    nR = 10  # Number of Radii
   

    # Azimuthal points
    nth = 15

   #given the radii, each circle has nth points 
   #points on the boundary
   
    # allocate space
    ldim = 3
    npoints = len(nz)*nth
    
    print('Allocated {0} points'.format(npoints))



pwal = [0] * npoints

k=0
#for testing
for il in range(len(nz)):
     for im in reversed(range(nth)):
         pwal[k]=(il+1)*nR*nth-im
         k=k+1


    #for testing
if __name__ == "__main__":
    fname = 'int_fld'
    data = read_int_fld(fname)
fptr_MD_2 = open("WSSdata.txt", "w")
for il in range(len(pwal)-1):

     lptn = data.pset[pwal[il]]
     data_pos = getattr(lptn,'pos')
        
     data_sts = getattr(lptn,'stat')
        
     mu=1.0
     rl=np.sqrt(data_pos[0]**2+data_pos[1]**2)
     if rl==0:
        nx=data_pos[0]
        ny=data_pos[1]
     else:
        nx=data_pos[0]/rl
        ny=data_pos[1]/rl
     tx=0.5*mu*((data_sts[78]+data_sts[78])*nx+(data_sts[79]+data_sts[81])*ny)
     ty=0.5*mu*((data_sts[82]+data_sts[82])*ny+(data_sts[79]+data_sts[81])*nx)
     tz=0.5*mu*((data_sts[80]+data_sts[84])*nx+(data_sts[83]+data_sts[85])*ny)
     data_T=[data_pos, tx, ty]
        
     datastr = "{} {} {} {} {} {}\n".format(float(data_pos[0]), float((data_pos[1])),float(data_pos[2]), float((tx)), float((ty),float((tz)))
     
     fptr_MD_2.writelines(datastr)
fptr_MD_2.close()




