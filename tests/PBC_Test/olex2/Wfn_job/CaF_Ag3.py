#!/usr/bin/env python

import numpy

TYPE_MAP = [
  [1],  # S
  [2, 3, 4],  # P
  [5, 8, 9, 6, 10, 7],  # D
  [11,14,15,17,20,18,12,16,19,13],  # F
  [21,24,25,30,33,31,26,34,35,28,22,27,32,29,23],  # G
  [56,55,54,53,52,51,50,49,48,47,46,45,44,43,42,41,40,39,38,37,36],  # H
]

def write_wfn(fout, mol, mo_coeff, mo_energy, mo_occ, tot_ener):
  from pyscf.x2c import x2c
  mol, ctr = x2c._uncontract_mol(mol, True, 0.)
  mo_coeff = numpy.dot(ctr, mo_coeff)

  nmo = mo_coeff.shape[1]
  mo_cart = []
  from pyscf import gto
  centers = []
  types = []
  exps = []
  p0 = 0
  for ib in range(mol.nbas):
    ia = mol.bas_atom(ib)
    l = mol.bas_angular(ib)
    es = mol.bas_exp(ib)
    c = mol._libcint_ctr_coeff(ib)
    np, nc = c.shape
    nd = nc*(2*l+1)
    mosub = mo_coeff[p0:p0+nd].reshape(-1,nc,nmo)
    c2s = gto.cart2sph(l)
    mosub = numpy.einsum('yki,cy,pk->pci', mosub, c2s, c)
    mo_cart.append(mosub.transpose(1,0,2).reshape(-1,nmo))

    for t in TYPE_MAP[l]:
        types.append([t]*np)
    ncart = mol.bas_len_cart(ib)
    exps.extend([es]*ncart)
    centers.extend([ia+1]*(np*ncart))
    p0 += nd
  mo_cart = numpy.vstack(mo_cart)
  centers = numpy.hstack(centers)
  types = numpy.hstack(types)
  exps = numpy.hstack(exps)
  nprim, nmo = mo_cart.shape

  fout.write('From PySCF\\n')
  fout.write('GAUSSIAN %14d MOL ORBITALS %6d PRIMITIVES %8d NUCLEI\\n'%(mo_cart.shape[1], mo_cart.shape[0], mol.natm))
  for ia in range(mol.natm):
    x, y, z = mol.atom_coord(ia)
    fout.write('%3s%8d (CENTRE%3d) %12.8f%12.8f%12.8f  CHARGE = %4.1f\\n'%(mol.atom_pure_symbol(ia), ia+1, ia+1, x, y, z, mol.atom_charge(ia)))
  for i0, i1 in lib.prange(0, nprim, 20):
    fout.write('CENTRE ASSIGNMENTS  %s\\n'% ''.join('%3d'%x for x in centers[i0:i1]))
  for i0, i1 in lib.prange(0, nprim, 20):
    fout.write('TYPE ASSIGNMENTS    %s\\n'% ''.join('%3d'%x for x in types[i0:i1]))
  for i0, i1 in lib.prange(0, nprim, 5):
    fout.write('EXPONENTS  %s\\n'% ' '.join('%13.7E'%x for x in exps[i0:i1]))

  for k in range(nmo):
      mo = mo_cart[:,k]
      fout.write('MO  %-12d          OCC NO = %12.8f ORB. ENERGY = %12.8f\\n'%(k+1, mo_occ[k], mo_energy[k]))
      for i0, i1 in lib.prange(0, nprim, 5):
        fout.write(' %s\\n' % ' '.join('%15.8E'%x for x in mo[i0:i1]))
  fout.write('END DATA\\n')
  fout.write(' THE SCF ENERGY =%20.12f THE VIRIAL(-V/T)=   0.00000000\\n'%tot_ener)

from pyscf.pbc import gto, scf, dft, df
from pyscf import lib
lib.num_threads(6)
cell = gto.M(
  atom = '''Ca  2.7254760000 2.7254760000 2.7254760000
Ca  2.7254760000 0.0000000000 0.0000000000
Ca  0.0000000000 2.7254760000 0.0000000000
F  4.0882140000 4.0882140000 4.0882140000
F  1.3627380000 4.0882140000 1.3627380000
F  4.0882140000 1.3627380000 1.3627380000
F  1.3627380000 1.3627380000 4.0882140000
F  4.0882140000 4.0882140000 1.3627380000
F  4.0882140000 1.3627380000 4.0882140000
F  1.3627380000 4.0882140000 4.0882140000
Ca  0.0000000000 0.0000000000 2.7254760000
F  1.3627380000 1.3627380000 1.3627380000
''',
  verbose = 5,
)
cell.output = 'CaF_Ag3_pyscf.log'
cell.a = '''5.45095 0.00000 0.00000
0.00000 5.45095 0.00000
0.00000 0.00000 5.45095'''
cell.charge = 0
cell.spin = 0
cell.max_memory = 15360.0
cell.precision = 1.0e-06
cell.basis = {'Ca': [[0,
                (1915.4348000, 0.0646240),
                (289.5332400, 0.3798380),
                (63.1063520, 0.6783290),],
[0,
                (80.3974400, -0.1093030),
                (17.3307500, 0.1089000),
                (5.0836240, 0.9492770),],
[0,
                (4.7822290, -0.2816070),
                (1.4625580, 0.3410510),
                (0.4792230, 0.8381040),],
[0,
                (0.4396820, -0.2697050),
                (0.0591300, 1.1132930),],
[0,
                (0.0238970, 1.0000000),],
[1,
                (80.3974400, 0.1354330),
                (17.3307500, 0.5372220),
                (5.0836240, 0.5018040),],
[1,
                (4.7822290, 0.0190090),
                (1.4625580, 0.4360380),
                (0.4792230, 0.6386710),],
[1,
                (0.4396820, 0.0003080),
                (0.0591300, 0.9998960),],
[1,
                (0.0238970, 1.0000000),],
],
'F': [[0,
                (413.8010000, 0.0585483),
                (62.2446000, 0.3493080),
                (13.4340000, 0.7096320),],
[0,
                (9.7775900, -0.4073270),
                (2.0861700, 1.2231400),],
[0,
                (0.4823830, 1.0000000),],
[1,
                (9.7775900, 0.2466800),
                (2.0861700, 0.8523210),],
[1,
                (0.4823830, 1.0000000),],
],

}
cell.exp_to_discard = 0.1
cell.build()
cf = dft.RKS(cell)
cf.xc = 'pbe,pbe'
cf = cf.mix_density_fit()
cf.grids.level = 0

#gdf = df.GDF(cell)
#cf.with_df = gdf
cf.with_df.auxbasis = 'def2-universal-jkfit'
cf.diis_space = 10
cf.conv_tol = 0.0033
cf.conv_tol_grad = 1e-2
cf.level_shift = 0.25
cf.damp = 0.7
cf.kernel()
print("Switching to SOSCF and shutting down damping & levelshift")
cf.conv_tol = 1e-9
cf.conv_tol_grad = 1e-5
cf.level_shift = 0.0
cf.damp = 0.0
cf = scf.newton(cf)
ener = cf.kernel()
with open('CaF_Ag3.wfn', 'w') as f1:
  write_wfn(f1,cell,cf.mo_coeff,cf.mo_energy,cf.mo_occ,ener)