# Example script for reading int_fld data
#    data.pord = pord
    data.start_time = stime[0]
    data.end_time = stime[1]
    data.effav_time = stime[2]
    data.dt = stime[3]
    data.nrec = nrec[0]
    data.int_time = itime[0]

    # read coordinates
    for il in range(data.npoints)
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

    # read derivatives
    for il in range(data.nderiv):
        for jl in range(data.npoints):
            ldrv = read_flt(infile,emode,wdsize,1)
            lptn = data.pset[jl]
            data_drv = getattr(lptn,'deriv')
            data_drv[il] =  ldrv[0]
    return data

def print_sim_data(data):
    """print simulation data"""
    print('Simulation data:')
    print('    Re = {0}'.format(data.re))
    print('    Lx = {0}'.format(data.bsize[0]))
    print('    Ly = {0}'.format(data.bsize[1]))
    print('    Lz = {0}'.format(data.bsize[2]))
    print('    nelx = {0}'.format(data.belem[0]))
    print('    nely = {0}'.format(data.belem[1]))
    print('    nelz = {0}'.format(data.belem[2]))
    print('    lx1 = {0}'.format(data.pord[0]))
    print('    ly1 = {0}'.format(data.pord[1]))
    print('    lz1 = {0}'.format(data.pord[2]))
    print('    nstat = {0}'.format(data.nstat))
    print('    nderiv = {0}'.format(data.nderiv))
    print('    start time = {0}'.format(data.start_time))
    print('    end time = {0}'.format(data.end_time))
    print('    effective average time = {0}'.format(data.effav_time))
    print('    time step = {0}'.format(data.dt))
    print('    nrec = {0}'.format(data.nrec))
    print('    time interval = {0}'.format(data.int_time))
    print('    npoints = {0}'.format(data.npoints))

def print_point_data(data,il):
    """print data related to single point"""
    print('Point data, npt = {0}'.format(il+1))
    lptn = data.pset[il]
    data_pos = getattr(lptn,'pos')
    print(data_pos)
    data_sts = getattr(lptn,'stat')
    print(data_sts)
    data_drv = getattr(lptn,'deriv')
    print(data_drv.size)
if __name__ == "__main__":
    fname = 'int_fld'
    data = read_int_fld(fname)

    print_sim_data(data)



    # place for operations

