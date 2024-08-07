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

  fout.write('From PySCF\n')
  fout.write('GAUSSIAN %14d MOL ORBITALS %6d PRIMITIVES %8d NUCLEI\n'%(mo_cart.shape[1], mo_cart.shape[0], mol.natm))
  for ia in range(mol.natm):
    x, y, z = mol.atom_coord(ia)
    fout.write('%3s%8d (CENTRE%3d) %12.8f%12.8f%12.8f  CHARGE = %4.1f\n'%(mol.atom_pure_symbol(ia), ia+1, ia+1, x, y, z, mol.atom_charge(ia)))
  for i0, i1 in lib.prange(0, nprim, 20):
    fout.write('CENTRE ASSIGNMENTS  %s\n'% ''.join('%3d'%x for x in centers[i0:i1]))
  for i0, i1 in lib.prange(0, nprim, 20):
    fout.write('TYPE ASSIGNMENTS    %s\n'% ''.join('%3d'%x for x in types[i0:i1]))
  for i0, i1 in lib.prange(0, nprim, 5):
    fout.write('EXPONENTS  %s\n'% ' '.join('%13.7E'%x for x in exps[i0:i1]))

  for k in range(nmo):
      mo = mo_cart[:,k]
      fout.write('MO  %-12d          OCC NO = %12.8f ORB. ENERGY = %12.8f\n'%(k+1, mo_occ[k], mo_energy[k]))
      for i0, i1 in lib.prange(0, nprim, 5):
        fout.write(' %s\n' % ' '.join('%15.8E'%x for x in mo[i0:i1]))
  fout.write('END DATA\n')
  fout.write(' THE SCF ENERGY =%20.12f THE VIRIAL(-V/T)=   0.00000000\n'%tot_ener)

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
                (2.402654E+06, 9.310000E-06),
                (3.597898E+05, 7.239000E-05),
                (8.187809E+04, 3.805900E-04),
                (2.319089E+04, 1.604530E-03),
                (7.565212E+03, 5.807800E-03),
                (2.730702E+03, 1.859566E-02),
                (1.064640E+03, 5.287776E-02),
                (4.410605E+02, 1.301515E-01),
                (1.917269E+02, 2.593147E-01),
                (8.653774E+01, 3.614961E-01),
                (3.989924E+01, 2.641116E-01),
                (1.764065E+01, 5.709398E-02),
                (8.359990E+00, -1.822000E-03),
                (3.951330E+00, 2.111640E-03),
                (1.713400E+00, -9.770900E-04),
                (8.108600E-01, 4.558100E-04),
                (3.602500E-01, -1.914600E-04),
                (8.108000E-02, 9.171000E-05),
                (4.484000E-02, -7.874000E-05),
                (2.143000E-02, 2.232000E-05),],
[0,
                (2.402654E+06, -2.700000E-06),
                (3.597898E+05, -2.102000E-05),
                (8.187809E+04, -1.105200E-04),
                (2.319089E+04, -4.666600E-04),
                (7.565212E+03, -1.695130E-03),
                (2.730702E+03, -5.483460E-03),
                (1.064640E+03, -1.596597E-02),
                (4.410605E+02, -4.151384E-02),
                (1.917269E+02, -9.286387E-02),
                (8.653774E+01, -1.653166E-01),
                (3.989924E+01, -1.766407E-01),
                (1.764065E+01, 6.444422E-02),
                (8.359990E+00, 5.108796E-01),
                (3.951330E+00, 4.946380E-01),
                (1.713400E+00, 8.750094E-02),
                (8.108600E-01, -3.591010E-03),
                (3.602500E-01, 2.492200E-03),
                (8.108000E-02, -7.582500E-04),
                (4.484000E-02, 6.458700E-04),
                (2.143000E-02, -1.817400E-04),],
[0,
                (2.402654E+06, 9.300000E-07),
                (3.597898E+05, 7.250000E-06),
                (8.187809E+04, 3.811000E-05),
                (2.319089E+04, 1.610100E-04),
                (7.565212E+03, 5.846600E-04),
                (2.730702E+03, 1.895040E-03),
                (1.064640E+03, 5.525250E-03),
                (4.410605E+02, 1.447010E-02),
                (1.917269E+02, 3.271581E-02),
                (8.653774E+01, 6.003190E-02),
                (3.989924E+01, 6.701683E-02),
                (1.764065E+01, -2.593660E-02),
                (8.359990E+00, -2.674747E-01),
                (3.951330E+00, -4.269726E-01),
                (1.713400E+00, 6.796405E-02),
                (8.108600E-01, 7.102054E-01),
                (3.602500E-01, 4.418006E-01),
                (8.108000E-02, 2.193139E-02),
                (4.484000E-02, -1.186924E-02),
                (2.143000E-02, 2.652720E-03),],
[0,
                (2.402654E+06, -2.200000E-07),
                (3.597898E+05, -1.730000E-06),
                (8.187809E+04, -9.100000E-06),
                (2.319089E+04, -3.844000E-05),
                (7.565212E+03, -1.396500E-04),
                (2.730702E+03, -4.525100E-04),
                (1.064640E+03, -1.320380E-03),
                (4.410605E+02, -3.458410E-03),
                (1.917269E+02, -7.836860E-03),
                (8.653774E+01, -1.441517E-02),
                (3.989924E+01, -1.621648E-02),
                (1.764065E+01, 6.344200E-03),
                (8.359990E+00, 6.740626E-02),
                (3.951330E+00, 1.144740E-01),
                (1.713400E+00, -2.634477E-02),
                (8.108600E-01, -2.336989E-01),
                (3.602500E-01, -3.160754E-01),
                (8.108000E-02, 3.328194E-01),
                (4.484000E-02, 5.611103E-01),
                (2.143000E-02, 2.808178E-01),],
[0,
                (8.108000E-02, 1.000000E+00),],
[0,
                (2.143000E-02, 1.000000E+00),],
[1,
                (4.061289E+03, 1.979900E-04),
                (9.622465E+02, 1.732080E-03),
                (3.121686E+02, 9.533790E-03),
                (1.187144E+02, 3.839012E-02),
                (4.980670E+01, 1.167588E-01),
                (2.225998E+01, 2.562687E-01),
                (1.028764E+01, 3.797808E-01),
                (4.861154E+00, 3.082933E-01),
                (2.248773E+00, 8.592090E-02),
                (1.033662E+00, 2.120670E-03),
                (4.641320E-01, 1.288800E-03),
                (1.987500E-01, -4.683500E-04),
                (6.739000E-02, 1.472800E-04),
                (2.542000E-02, -5.288000E-05),],
[1,
                (4.061289E+03, -6.455000E-05),
                (9.622465E+02, -5.645800E-04),
                (3.121686E+02, -3.131250E-03),
                (1.187144E+02, -1.274086E-02),
                (4.980670E+01, -3.991403E-02),
                (2.225998E+01, -9.050448E-02),
                (1.028764E+01, -1.426190E-01),
                (4.861154E+00, -1.098090E-01),
                (2.248773E+00, 1.516249E-01),
                (1.033662E+00, 4.617641E-01),
                (4.641320E-01, 4.326003E-01),
                (1.987500E-01, 1.112741E-01),
                (6.739000E-02, 2.528740E-03),
                (2.542000E-02, 7.103200E-04),],
[1,
                (4.061289E+03, 1.336000E-05),
                (9.622465E+02, 1.185200E-04),
                (3.121686E+02, 6.487200E-04),
                (1.187144E+02, 2.679930E-03),
                (4.980670E+01, 8.285130E-03),
                (2.225998E+01, 1.921235E-02),
                (1.028764E+01, 2.954984E-02),
                (4.861154E+00, 2.431232E-02),
                (2.248773E+00, -4.111230E-02),
                (1.033662E+00, -1.041976E-01),
                (4.641320E-01, -1.503654E-01),
                (1.987500E-01, 2.431729E-02),
                (6.739000E-02, 5.986116E-01),
                (2.542000E-02, 4.813557E-01),],
[1,
                (1.987500E-01, 1.000000E+00),],
[1,
                (2.542000E-02, 1.000000E+00),],
[2,
                (1.694623E+01, 1.547300E-02),
                (4.472120E+00, 7.887400E-02),
                (1.438090E+00, 2.087800E-01),
                (4.669900E-01, 3.302130E-01),
                (1.415100E-01, 4.375620E-01),
                (4.164000E-02, 3.747900E-01),],
[2,
                (1.415100E-01, 1.000000E+00),],
[2,
                (4.164000E-02, 1.000000E+00),],
[3,
                (1.509000E-01, 1.0000000),],
],
'F': [[0,
                (19500.00000, 0.5070000000E-03),
                (2923.000000, 0.3923000000E-02),
                (664.5000000, 0.2020000000E-01),
                (187.5000000, 0.7901000000E-01),
                (60.62000000, 0.2304390000),
                (21.42000000, 0.4328720000),
                (7.950000000, 0.3499640000),
                (0.8815000000, -0.7892000000E-02),],
[0,
                (19500.00000, -0.1170000000E-03),
                (2923.000000, -0.9120000000E-03),
                (664.5000000, -0.4717000000E-02),
                (187.5000000, -0.1908600000E-01),
                (60.62000000, -0.5965500000E-01),
                (21.42000000, -0.1400100000),
                (7.950000000, -0.1767820000),
                (0.8815000000, 0.6050430000),],
[0,
                (2.257000000, 1.000000000),],
[0,
                (0.3041000000, 1.000000000),],
[1,
                (43.88000000, 0.1666500000E-01),
                (9.926000000, 0.1044720000),
                (2.930000000, 0.3172600000),],
[1,
                (0.9132000000, 1.000000000),],
[1,
                (0.2672000000, 1.000000000),],
[2,
                (3.107000000, 1.000000000),],
[2,
                (0.8550000000, 1.000000000),],
[3,
                (1.917000000, 1.000000000),],
],

}
cell.exp_to_discard = 0.1
cell.ke_cutoff=100
cell.build()
cf = dft.RKS(cell)
cf.xc = 'pbe,pbe'
#cf = cf.mix_density_fit()
cf.grids.level = 0
mdf = df.MDF(cell)
mdf.linear_dep_threshold = 1e-6
mdf.mesh = [14,14,14]
cf.with_df = mdf
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
with open('CaF_cc_pVTZ.wfn', 'w') as f1:
  write_wfn(f1,cell,cf.mo_coeff,cf.mo_energy,cf.mo_occ,ener)