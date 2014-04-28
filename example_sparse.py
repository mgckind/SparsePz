#!/usr/bin/env python
__author__ = 'Matias Carrasco Kind'
from numpy import *
import pdf_storage as ps
import sys, os
try:
    from mpi4py import MPI
    PLL = 'MPI'
except ImportError:
    PLL = 'SERIAL'
import pyfits as pf


#This function is borrowed from MLZ utils
def get_limits(ntot, Nproc, rank):
    """
    Get limits for farming an array to multiple processors

    :param int ntot: Number of objects in array
    :param int Nproc: number of processor
    :param int rank: current processor id
    :return: L1,L2 the limits of the array for given processor
    :rtype: int, int
    """
    jpproc = zeros(Nproc) + int(ntot / Nproc)
    for i in xrange(Nproc):
        if (i < ntot % Nproc): jpproc[i] += 1
    jpproc = map(int, jpproc)
    st = rank
    st = sum(jpproc[:rank]) - 1
    s0 = int(st + 1)
    s1 = int(st + jpproc[rank]) + 1
    return s0, s1

if PLL=='MPI':
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
else:
    size=1
    rank=0
    
filein = 'CFHTLens_sample.P.npy'
#FORMAT FILE, EACH ROW IS THE PDF FOR EACH GALAXY, LAST ROW IS THE REDSHIFT POSITION

P = load(filein)
Ntot = len(P)-1 #last row is redshift

if rank == 0:
    print "Total Galaxies = ", Ntot
z = P[-1]
dz = z[1] - z[0]
if rank == 0:
    print 'dz = ', dz

mu = [min(z), max(z)]
Nmu = 250 #len(z)
max_sig = (max(z) - min(z)) / 12.
min_sig = dz / 6.
Nsig = int(ceil(2. * (max_sig - min_sig) / dz))
sig = [min_sig, max_sig]
Nv = 3
Nsig = 80
NA = Nmu * Nsig * Nv


if rank == 0:
    print 'Nmu, Nsig, Nv = ', '[', Nmu, ',', Nsig, ',', Nv, ']'
    print 'Total bases in dictionary', NA

#Create dictionary
if rank == 0:
    print 'Creating Dictionary...'
A = ps.create_voigt_dict(z, mu, Nmu, sig, Nsig, Nv)
bigD = {}

toler = 1.e-10
Nsparse = 20
Ncoef = 32001
AA = linspace(0, 1, Ncoef)
Da = AA[1] - AA[0]

if rank==0:
    print 'Nsparse (number of bases) = ',Nsparse

bigD['z'] = z
bigD['mu'] = mu
bigD['sig'] = sig
bigD['dims'] = [Nmu, Nsig, Nv, Ncoef]
bigD['Nsparse'] = Nsparse
bigD['Ntot']=Ntot

if rank == 0:
    for i in xrange(size):
        Xs_0, Xs_1 = get_limits(Ntot, size, i)
        print Xs_0, ' ', Xs_1, ' -------------> to core ', i

s0, s1 = get_limits(Ntot, size, rank)
P = P[s0:s1]

if rank == 0:
    print 'Creating Sparse representation...'

for ik in xrange(Ntot):
    k = s0 + ik
    bigD[k] = {}
    try:
        pdf0 = P[ik]
    except:
        continue
    if sum(pdf0) > 0:
        pdf0 /= sum(pdf0)
    else:
        continue
    np=Nsparse
    Dind, Dval = ps.sparse_basis(A, pdf0, np)
    if len(Dind) <=1 : continue
    bigD[k]['sparse'] = [Dind, Dval]
    if max(Dval) > 0 :
        Dvalm = Dval / max(Dval)
        index = array(map(round, (Dvalm / Da)), dtype='int')
    else:
        index=zeros(len(Dind))
    bigD[k]['sparse_ind'] = array(map(ps.combine_int, index, Dind))

    #swap back columns
    A[:, [Dind]] = A[:, [arange(len(Dind))]]


print 'Done with processor: ',rank
if PLL=='MPI':comm.Barrier()

if PLL=='MPI':
    if rank == 0:
        for srank in xrange(1,size):
            temp=comm.recv(source=srank, tag=srank*2)
            bigD.update(temp)
            del temp
    else:
        comm.send(bigD,dest=0,tag=rank*2)
    comm.Barrier()
    
if rank==0:
    print 'Writing fits file (example_out.fits)'
    ALL=zeros((Ntot,Nsparse),dtype='int')
    for i in xrange(Ntot):
        if bigD.has_key(i):
            idd=bigD[i]['sparse_ind']
            ALL[i,0:len(idd)]=idd
    head=pf.Header()
    head['N_TOT']=Ntot
    head['N_MU']=bigD['dims'][0]
    head['N_SIGMA']=bigD['dims'][1]
    head['N_VOIGT']=bigD['dims'][2]
    head['N_COEF']=bigD['dims'][3]
    head['N_SPARSE']=bigD['Nsparse']
    head['MU1']=bigD['mu'][0]
    head['MU2']=bigD['mu'][1]
    head['SIGMA1']=bigD['sig'][0]
    head['SIGMA2']=bigD['sig'][1]
    col1=pf.Column(name='redshift',format='E',array=bigD['z'])
    fmt='%dJ' % bigD['Nsparse']
    col2=pf.Column(name='Sparse_indices',format=fmt,array=ALL)
    table1=pf.new_table(pf.ColDefs([col1]))
    table2=pf.new_table(pf.ColDefs([col2]))
    prihdu=pf.PrimaryHDU(header=head)
    hdulist = pf.HDUList([prihdu,table1,table2])
    hdulist.writeto('example_out.fits',clobber=True)

if PLL=='MPI': 
    comm.Barrier()
    MPI.Finalize()






