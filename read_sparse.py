#!/usr/bin/env python
__author__ = 'Matias Carrasco Kind'
from numpy import *
import pyfits as pf
import matplotlib.pyplot as plt
import pdf_storage as ps
import random as rn


PO=load('CFHTLens_sample.P.npy')

F=pf.open('example_out.fits')
H=F[0].header
Ntot=H['N_TOT']
Nmu=H['N_MU']
Nsig=H['N_SIGMA']
Nv=H['N_VOIGT']
Ncoef=H['N_COEF']
Nspa=H['N_SPARSE']
mu=[H['MU1'],H['MU2']]
sig=[H['SIGMA1'],H['SIGMA2']]
z=F[1].data.field('redshift')
P=F[2].data.field('Sparse_indices')

VALS=linspace(0,1,Ncoef)
dVals=VALS[1]-VALS[0]

k=rn.sample(xrange(len(PO)-1),1)[0]


sp_ind=array(map(ps.get_N,P[k]))

spi=sp_ind[:,0]
Dind2=sp_ind[:,1]
vals=spi*dVals
rep_pdf = ps.reconstruct_pdf_v(Dind2, vals, z, mu, Nmu, sig, Nsig, Nv)

plt.plot(z,PO[k]/sum(PO[k]),label='original')
plt.plot(z,rep_pdf,label='Sparse rep')
plt.xlabel('redshift')
plt.ylabel('P(z)')
plt.legend(loc=0)
title='Galaxy example No: %d' % k
plt.title(title)

plt.show()

