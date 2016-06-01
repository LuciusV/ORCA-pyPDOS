#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Partial density-of-states (PDOS) plotter for ORCA (https://orcaforum.cec.mpg.de/) written in python
AUTHOR: Dr. Evgeny V. Tikhonov, Moscow State University, Physics department, email: e.tikhonov@physics.msu.ru
"""

"""
N-th molecular orbital represented as:
\Psi(N) = \sum_m C_{mn}\phi_m
where \phi_m is m-th atomic orbital
It is sliced as Cmn[N,:]
"""

"""
TODO (in order of importance:
1. Interactive and not-interactive modes
2. Adding __doc__ function, optimizing import from pylab 


KNOWN BUGS

(#001) For lines where there is no space between columns, like:
  3H   2s        -3.975631  0.012837 -2.950634 -0.009080-16.567733 -6.648447
parser fails to split line correctly (fixed)
MAY HAPPEN IF NO - SIGN WOULD BE BETWEEN TWO COLUMNS WITHOUT WHITESPACE
"""


from pylab import *
from scipy.stats import norm
from scipy.constants import codata
import argparse
import sys
import re

parser = argparse.ArgumentParser(description='Parse ORCA out file and plot partial DOS')
group = parser.add_mutually_exclusive_group()
parser.add_argument('filename', type=str, help="path to input file(s)", nargs='+')
parser.add_argument('-s', '--smear', type=float, default=0.1, help="smearing width (in eV)")
parser.add_argument('-d', '--dump', type=bool, default=False, action='store_true', help='store data in txt files')
parser.add_argument('-v', '--verbosity', action='count', default=0, help="verbosity level(-v,-vv,-vvv)")
group.add_argument('-u', '--unique', default=False, action='store_true', help="plot DOS for unique atoms")
group.add_argument('-a', '--atoms', type=int, help="array of atoms for which to plot DOS", nargs='+')
parser.add_argument('-L', '--lowest', type=float, default=-20.0, help="Lowest energy range for plot (in eV)")
parser.add_argument('-H', '--highest', type=float, default=0.0, help="Highest energy range for plot (in eV)")
args = parser.parse_args()

atoms_listed = False

if args.atoms:
    args.atoms = unique(array(args.atoms))
    if args.verbosity >= 3: print type(args.atoms), size(args.atoms), args.atoms
    atoms_listed = True

def gaussian(x, mu, sig):
    return exp(-(x - mu)*(x - mu)/(2*sig*sig))

def list_to_array(List):
    arr = []
    for line in List:
        for x in line[:-1].split():
            arr.append(float(x))
    return array(arr)

def find_Ef(eigs,ocns):
    global HFTyp
    if HFTyp == 'RHF':
        e = eigs[0]
        f = ocns[0]
        for ei,fi in zip(e,f):
            if fi < 0.5:
                LUMO = ei
                break
            HOMO = ei
        return (HOMO+LUMO)*0.5
    if HFTyp == 'UHF' :
        e_a = eigs[0]
        e_b = eigs[1]
        f_a = ocns[0]
        f_b = ocns[1]
        for ei,fi in zip(e_a,f_a):
            if fi < 0.5:
                LUMO_A = ei
                break
            HOMO_A = ei
        for ei,fi in zip(e_b,f_b):
            if fi < 0.5:
                LUMO_B = ei
                break
            HOMO_B = ei
        HOMO = max(HOMO_A,HOMO_B)
        LUMO = min(LUMO_A,LUMO_B)
        return (HOMO+LUMO)*0.5

def get_mo_decomposition(mo_number, data, renorm = False):
    global DIM
    global HFTyp
    ln = 1
    c = []
    index = mo_number%6 + 2
    shift = (DIM+4)*(mo_number/6)
    for line in data:
        if ln >= 5 + shift and ln <= (DIM + 4) + shift:
            tmp = float(line[:-1].split()[index])
            c.append(tmp)
        ln += 1
    arr = array(c)
    if renorm and HFTyp == 'RHF':
        return 2*arr*arr/(arr*arr).sum()
    if renorm and HFTyp == 'UHF':
        return arr*arr/(arr*arr).sum()
    else:
        return arr

def main():
    global DIM
    global HFTyp
    for f in args.filename:
        if args.verbosity >= 1: print "Loading file "+f
        species_found = False
        species_amount_found = False
        hftyp_found = False
        nel_found = False
        dim_found = False
        mol_offset_found = False
        if args.verbosity >= 2: print "Plot range: ["+str(args.lowest)+" "+str(args.highest)+"] eV"
        with open(f,'r') as log:
            ln = 1
            data = []
            data_a = []
            data_b = []
            species = []
            for line in log:
                if not species_found:
                    if "BASIS SET INFORMATION" in line:
                        species_offset = ln
                        species_found = True
                if species_found and ln == species_offset + 2:
                    species_amount = int(line[:-1].split()[2])
                    species_amount_found = True
                if species_amount_found and ln >= species_offset + 4 and ln <= species_offset + 4 + species_amount - 1:
                    species.append(line[:-1].split()[3])
                if not hftyp_found:
                    if "Hartree-Fock type" in line:
                        HFTyp = str(line.split()[-1])
                        if HFTyp == 'UHF':
                            del(data)
                        if HFTyp == 'RHF':
                            del(data_a)
                            del(data_b)
                        hftyp_found = True
                if not nel_found:
                    if "Number of Electrons" in line:
                        NEL = int(line.split()[-1])
                        nel_found = True
                if not dim_found:
                    if "Basis Dimension" in line:
                        DIM = int(line.split()[-1])
                        dim_found = True
                        if DIM%6 != 0: shift = 1
                        else: shift = 0
                if not mol_offset_found:
                    if "MOLECULAR ORBITALS" in line:
                        MOL_OFFSET = ln
                        mol_offset_found = True
                if mol_offset_found and HFTyp == 'RHF' and ln >= MOL_OFFSET+2 and ln <= MOL_OFFSET+1+(DIM+4)*(DIM/6+shift):
                    expr = re.search("[0-9]-[0-9]",line) # to fix bug #001
                    if expr:
                        match = expr.group(0)[0]+" -"+expr.group(0)[-1]
                        newline = re.sub("[0-9]-[0-9]", match, line)
                        data.append(newline[:-1])
                    else:
                        data.append(line[:-1])
                if mol_offset_found and HFTyp == 'UHF' and ln >= MOL_OFFSET+2 and ln <= MOL_OFFSET+1+(DIM+4)*(DIM/6+shift):
                    expr = re.search("[0-9]-[0-9]",line)
                    if expr:
                        match = expr.group(0)[0]+" -"+expr.group(0)[-1]
                        newline = re.sub("[0-9]-[0-9]", match, line)
                        data_a.append(newline[:-1])
                    else:
                        data_a.append(line[:-1])
                if mol_offset_found and HFTyp == 'UHF' and ln >= MOL_OFFSET+2+(DIM+4)*(DIM/6+shift) +1 and ln <= MOL_OFFSET+1+2*(DIM+4)*(DIM/6+shift) +1:
                    expr = re.search("[0-9]-[0-9]",line)
                    if expr:
                        match = expr.group(0)[0]+" -"+expr.group(0)[-1]
                        newline = re.sub("[0-9]-[0-9]", match, line)
                        data_b.append(newline[:-1])
                    else:
                        data_b.append(line[:-1])
                ln += 1
            if not mol_offset_found:
                sys.exit("Molecular orbitals not found. Please rerun ORCA job with option \"%output print [p_mos] 1 end\".")
            BASISAO = []
            if atoms_listed:
                ATOMLIST = []
            if HFTyp == 'RHF':
                for line in data[4:4+DIM]:
                    if args.unique:
                        BASISAO.append(line.split()[0])
                    elif atoms_listed:
                        AO = line.split()[0]
                        BASISAO.append(AO)
                        anumber = re.sub("[a-zA-Z]","",AO)
                        if intersect1d(array(anumber), args.atoms).size:
                            ATOMLIST.append(AO)
                    else:
                        BASISAO.append(re.sub("[0-9]", "", line.split()[0]))
                if args.unique:
                    species = set(BASISAO)
                    species = list(species)
                    species.sort()
                if atoms_listed:
                    species = set(ATOMLIST)
                    species = list(species)
                    species.sort()
                if args.verbosity >= 2: print BASISAO
                if args.verbosity >= 1: print "Unique species:", species
                MOS_EIG = list_to_array(data[1::DIM+4])
                MOS_EIG = MOS_EIG*codata.value('Hartree energy in eV')
                MOS_OCC = list_to_array(data[2::DIM+4])
                fermilevel = find_Ef([MOS_EIG],[MOS_OCC])
                domain = linspace(args.lowest, args.highest, 10000)
                if atoms_listed:
                    full_occ_sum = zeros(domain.shape)
                    full_free_sum = zeros(domain.shape)
                Cmn = []
                for n in arange(DIM):
                    mo = get_mo_decomposition(n, data, True)
                    Cmn.append(mo)
                Cmn = array(Cmn)
                BASISAO = array(BASISAO)
                for atom in species:
                    Sum_occ = zeros(domain.shape)
                    Sum_free = zeros(domain.shape)
                    indices = where(BASISAO == atom)
                    for N,e in zip(arange(DIM),MOS_EIG):
                        if args.verbosity >= 3: print e,Cmn[N, indices].sum()
                        if MOS_OCC[N] > 0.5:
                            Sum_occ += Cmn[N, indices].sum()*gaussian(domain, e, args.smear)
                        else:
                            Sum_free += Cmn[N, indices].sum()*gaussian(domain, e, args.smear)
                    if atoms_listed:
                        full_occ_sum += Sum_occ 
                        full_free_sum += Sum_free
                    plot(domain, Sum_occ, label=atom+' occ.')
                    plot(domain, Sum_free,'--', label=atom+' free')
                    if args.dump:
                        np.savetxt(atom+'_occ.txt',Sum_occ)
                        np.savetxt(atom+'_free.txt',Sum_free)
                Sum_occ = zeros(domain.shape)
                Sum_free = zeros(domain.shape)
                for e,f in zip(MOS_EIG,MOS_OCC):
                    if f > 0.5:
                        Sum_occ += 2*gaussian(domain, e, args.smear)
                    else:
                        Sum_free += 2*gaussian(domain, e, args.smear)
                if atoms_listed:
                    plot(domain, full_occ_sum, color='#f1595f', lw = 2.0, label='full occ. of sel. atoms')
                    plot(domain, full_free_sum, '--', color='#f1595f', lw = 2.0, label='full occ. of sel. atoms')
                plot(domain, Sum_occ,'k', lw = 2.0, label='total occ. DOS')
                plot(domain, Sum_free,'k--', lw = 2.0, label='total free DOS')
                if args.dump:
                    np.savetxt('total_occ', Sum_occ)
                    np.savetxt('total_free', Sum_free)
            if HFTyp == 'UHF':
                for line in data_a[4:4+DIM]:
                    if args.unique:
                        BASISAO.append(line.split()[0])
                    elif atoms_listed:
                        AO = line.split()[0]
                        BASISAO.append(AO)
                        anumber = re.sub("[a-zA-Z]","",AO)
                        if intersect1d(array(anumber), args.atoms).size:
                            ATOMLIST.append(AO)
                    else:
                        BASISAO.append(re.sub("[0-9]", "", line.split()[0]))
                if args.unique:
                    species = set(BASISAO)
                    species = list(species)
                    species.sort()
                if atoms_listed:
                    species = set(ATOMLIST)
                    species = list(species)
                    species.sort()
                if args.verbosity >= 2: print BASISAO
                if args.verbosity >= 1: print "Unique species:", species
                MOS_EIG_A = list_to_array(data_a[1::DIM+4])
                MOS_EIG_A = MOS_EIG_A*codata.value('Hartree energy in eV')
                MOS_EIG_B = list_to_array(data_b[1::DIM+4])
                MOS_EIG_B = MOS_EIG_B*codata.value('Hartree energy in eV')
                MOS_OCC_A = list_to_array(data_a[2::DIM+4])
                MOS_OCC_B = list_to_array(data_b[2::DIM+4])
                fermilevel = find_Ef([MOS_EIG_A,MOS_EIG_B],[MOS_OCC_A,MOS_OCC_B])
                domain = linspace(args.lowest, args.highest, 10000)
                if atoms_listed:
                    full_occ_sum_a = zeros(domain.shape)
                    full_free_sum_a = zeros(domain.shape)
                    full_occ_sum_b = zeros(domain.shape)
                    full_free_sum_b = zeros(domain.shape)
                Cmn_a = []
                Cmn_b = []
                for n in arange(DIM):
                    mo_a = get_mo_decomposition(n, data_a, True)
                    mo_b = get_mo_decomposition(n, data_b, True)
                    Cmn_a.append(mo_a)
                    Cmn_b.append(mo_b)
                Cmn_a = array(Cmn_a)
                Cmn_b = array(Cmn_b)
                BASISAO = array(BASISAO)
                for atom in species:
                    Sum_a_occ = zeros(domain.shape)
                    Sum_a_free = zeros(domain.shape)
                    Sum_b_occ = zeros(domain.shape)
                    Sum_b_free = zeros(domain.shape)
                    indices = where(BASISAO == atom)
                    for N,e in zip(arange(DIM),MOS_EIG_A):
                        if args.verbosity >= 3: print e,Cmn_a[N, indices].sum()
                        if MOS_OCC_A[N] > 0.5:
                            Sum_a_occ += Cmn_a[N, indices].sum()*gaussian(domain, e, args.smear)
                        else:
                            Sum_a_free += Cmn_a[N, indices].sum()*gaussian(domain, e, args.smear)
                    for N,e in zip(arange(DIM),MOS_EIG_B):
                        if args.verbosity >=3 : print e,Cmn_b[N, indices].sum()
                        if MOS_OCC_B[N] > 0.5:
                            Sum_b_occ += Cmn_b[N, indices].sum()*gaussian(domain, e, args.smear)
                        else:
                            Sum_b_free += Cmn_b[N, indices].sum()*gaussian(domain, e, args.smear)
                    if atoms_listed:
                        full_occ_sum_a += Sum_a_occ
                        full_occ_sum_b += Sum_b_occ 
                        full_free_sum_a += Sum_a_free
                        full_free_sum_b += Sum_b_free
                    plot(domain, Sum_a_occ, label=atom+' up occ.')
                    plot(domain, Sum_a_free,'--', label=atom+' up free')
                    plot(domain, -Sum_b_occ, label=atom+' down occ.')
                    plot(domain, -Sum_b_free,'--', label=atom+' down free')
                    if args.dump:
                        np.savetxt(atom + '_up_occ.txt',Sum_a_occ)
                        np.savetxt(atom + '_up_free',Sum_a_free)
                        np.savetxt(atom + '_down_occ',Sum_b_occ)
                        np.savetxt(atom + '_down_free',Sum_b_free)
                Sum_a_occ = zeros(domain.shape)
                Sum_a_free = zeros(domain.shape)
                Sum_b_occ = zeros(domain.shape)
                Sum_b_free = zeros(domain.shape)
                for e,f in zip(MOS_EIG_A,MOS_OCC_A):
                    if f > 0.5:
                        Sum_a_occ += gaussian(domain, e, args.smear)
                    else:
                        Sum_a_free += gaussian(domain, e, args.smear)
                for e,f in zip(MOS_EIG_B,MOS_OCC_B):
                    if f > 0.5:
                        Sum_b_occ += gaussian(domain, e, args.smear)
                    else:
                        Sum_b_free += gaussian(domain, e, args.smear)
                if atoms_listed:
                    plot(domain, full_occ_sum_a, color='#f1595f', lw = 2.0, label='full occ. of sel. atoms')
                    plot(domain, full_free_sum_a, '--',color='#f1595f', lw = 2.0, label='full occ. of sel. atoms')
                    plot(domain, -full_occ_sum_b, color='#599ad3', lw = 2.0, label='full occ. of sel. atoms')
                    plot(domain, -full_free_sum_b, '--', color='#599ad3', lw = 2.0, label='full occ. of sel. atoms')
                plot(domain, Sum_a_occ, 'k', lw = 2.0, label='total up DOS')
                plot(domain, Sum_a_free, 'k--', lw = 2.0, label='total up DOS')
                plot(domain, -Sum_b_occ, 'k', lw = 2.0, label='total down DOS')
                plot(domain, -Sum_b_free, 'k--', lw = 2.0, label='total down DOS')
                if args.dump:
                    np.savetxt('total_up_occ', Sum_a_occ)
                    np.savetxt('total_up_free', Sum_a_free)
                    np.savetxt('total_down_occ', Sum_b_occ)
                    np.savetxt('total_down_free', Sum_b_free)
            axvline(x = fermilevel, ymin = 0.0, ymax = 1.0, color = 'k')
            legend(loc=9,ncol=6)
            show()
            log.close()

if __name__ == "__main__":
    main()
