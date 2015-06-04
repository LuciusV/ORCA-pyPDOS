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
1. Option of plotting PDOS for given array of atoms.
2. Interactive and not-interactive modes
3. Adding __doc__ function, optimizing import from pylab 

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
parser.add_argument('-v', '--verbosity', action='count', default=0, help="verbosity level(-v,-vv,-vvv)")
group.add_argument('-u', '--unique', default=False, action='store_true', help="plot DOS for unique atoms")
group.add_argument('-a', '--atoms', type=int, help="array of atoms for which to plot DOS", nargs='+')
parser.add_argument('-L', '--lowest', type=float, default=-20.0, help="Lowest energy range for plot (in eV)")
parser.add_argument('-H', '--highest', type=float, default=0.0, help="Highest energy range for plot (in eV)")
args = parser.parse_args()

atoms_listed = False

if args.atoms:
    args.atoms = array(args.atoms)
    atoms_listed = True

def gaussian(x, mu, sig):
    return exp(-(x - mu)*(x - mu)/(2*sig*sig))

def list_to_array(List):
    arr = []
    for line in List:
        for x in line[:-1].split():
            arr.append(float(x))
    return array(arr)

def plot_with_weight(value, weight):
    gauss = norm(loc = value, scale = args.smear)
    domain = linspace(value-args.smear*10, value+args.smear*10)
    plot(domain, gauss.pdf(domain))

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
                    data.append(line[:-1])
                if mol_offset_found and HFTyp == 'UHF' and ln >= MOL_OFFSET+2 and ln <= MOL_OFFSET+1+(DIM+4)*(DIM/6+shift):
                    data_a.append(line[:-1])
                if mol_offset_found and HFTyp == 'UHF' and ln >= MOL_OFFSET+2+(DIM+4)*(DIM/6+shift) +1 and ln <= MOL_OFFSET+1+2*(DIM+4)*(DIM/6+shift) +1:
                    data_b.append(line[:-1])
                ln += 1
            if not mol_offset_found:
                sys.exit("Molecular orbitals not found. Please rerun ORCA job with option \"%output print [p_mos] 1 end\".")
            BASISAO = []
            if HFTyp == 'RHF':
                for line in data[4:4+DIM]:
                    if args.unique:
                        BASISAO.append(line.split()[0])
                        species = set(BASISAO)
                        species = list(species)
                        species.sort()
                    elif atoms_listed:
                        BASISAO.append(re.sub("[0-9]", "", line.split()[0]))
                    else:
                        BASISAO.append(re.sub("[0-9]", "", line.split()[0]))
                if args.verbosity >= 2: print BASISAO
                if args.verbosity >= 1: print "Unique species:", species
                MOS_EIG = list_to_array(data[1::DIM+4])
                MOS_EIG = MOS_EIG*codata.value('Hartree energy in eV')
                MOS_OCC = list_to_array(data[2::DIM+4])
                domain = linspace(args.lowest, args.highest, 10000)
                Cmn = []
                for n in arange(DIM):
                    mo = get_mo_decomposition(n, data, True)
                    Cmn.append(mo)
                Cmn = array(Cmn)
                BASISAO = array(BASISAO)
                for atom in species:
                    Sum = zeros(domain.shape)
                    indices = where(BASISAO == atom)
                    for N,e in zip(arange(DIM),MOS_EIG):
                        if args.verbosity >= 3: print e,Cmn[N, indices].sum()
                        Sum += Cmn[N, indices].sum()*gaussian(domain, e, args.smear)
                    plot(domain, Sum, label=atom)
                Sum = zeros(domain.shape)
                for e in MOS_EIG:
                    Sum += 2*gaussian(domain, e, args.smear)
                plot(domain, Sum,'k',label='total DOS')
            if HFTyp == 'UHF':
                for line in data_a[4:4+DIM]:
                    if not args.unique:
                        BASISAO.append(re.sub("[0-9]", "", line.split()[0]))
                    else:
                        BASISAO.append(line.split()[0])
                        species = set(BASISAO)
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
                domain = linspace(args.lowest, args.highest, 10000)
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
                    Sum_a = zeros(domain.shape)
                    Sum_b = zeros(domain.shape)
                    indices = where(BASISAO == atom)
                    for N,e in zip(arange(DIM),MOS_EIG_A):
                        if args.verbosity >= 3: print e,Cmn_a[N, indices].sum()
                        Sum_a += Cmn_a[N, indices].sum()*gaussian(domain, e, args.smear)
                    for N,e in zip(arange(DIM),MOS_EIG_B):
                        if args.verbosity >=3 : print e,Cmn_b[N, indices].sum()
                        Sum_b += Cmn_b[N, indices].sum()*gaussian(domain, e, args.smear)
                    plot(domain, Sum_a, label=atom)
                    plot(domain, -Sum_b, label=atom)
                Sum_a = zeros(domain.shape)
                Sum_b = zeros(domain.shape)
                for e in MOS_EIG_A:
                    Sum_a += gaussian(domain, e, args.smear)
                for e in MOS_EIG_B:
                    Sum_b += gaussian(domain, e, args.smear)
                plot(domain, Sum_a, 'k', label='total DOS')
                plot(domain, -Sum_b, 'k', label='total DOS')
            legend()
            show()
            log.close()

if __name__ == "__main__":
    main()
