#!/usr/bin/env python
__author__ = 'Matias Carrasco Kind'
from numpy import *
import pyfits as pf
import matplotlib.pyplot as plt
import pdf_storage as ps
import random as rn


PO = load('CFHTLens_sample.P.npy')

F = pf.open('example_out.fits')
P = F[2].data.field('Sparse_indices')
F.close()

k = rn.sample(xrange(len(PO) - 1), 1)[0]
head = ps.read_header('example_out.fits')
z = head['z']

rep_pdf = ps.reconstruct_pdf_int(P[k], head)

plt.plot(z, PO[k] / sum(PO[k]), label='original')
plt.plot(z, rep_pdf, label='Sparse rep')
plt.xlabel('redshift')
plt.ylabel('P(z)')
plt.legend(loc=0)
title = 'Galaxy example No: %d' % k
plt.title(title)

plt.show()

