from ReadOrcaLog import ReadOrcaLog
from ReadMOComposition import ReadMOComposition
from PlotTotalDOS import PlotTotalDOS
from PlotPartialDOS import PlotPartialDOS
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import re

if __name__ == '__main__':

    parser = ArgumentParser(description='Parse ORCA out file and plot partial DOS')
    parser.add_argument('filename', type=str, help="path to input file(s)")
    parser.add_argument('-s', '--smear', type=float, default=0.1, help="smearing width (in eV)")
    parser.add_argument('-d', '--dump', default=False, action='store_true', help='store data in txt files')
    parser.add_argument('-i', '--interactive', default=False, action='store_true', help='display plots interactively')
    parser.add_argument('-v', '--verbosity', action='count', default=0, help="verbosity level(-v,-vv,-vvv)")
    parser.add_argument('-u', '--unique', default=False, action='store_true', help="also plot DOS for unique atoms")
    parser.add_argument('-a', '--atoms', type=int, help="array of atoms for which to plot DOS", nargs='+')
    parser.add_argument('-L', '--lowest', type=float, default=-20.0, help="Lowest energy range for plot (in eV)")
    parser.add_argument('-H', '--highest', type=float, default=0.0, help="Highest energy range for plot (in eV)")
    args = parser.parse_args()

    OrcaLog = args.filename

    f = open(OrcaLog, 'rb')
    Type, NumberOfElectrons, BasisDimension, Content, LenContent = ReadOrcaLog(f)
    f.close()

    if Type == 'RKS':
        Species, ContributionMatrix, Eigenvalues, Occupations = ReadMOComposition(Content, LenContent, Type,
                                                                                  BasisDimension)
        EigenvaluesB = None
        OccupationsB = None
        ContributionMatrixB = None

    elif Type == 'UKS':
        Species, ContributionMatrix, Eigenvalues, Occupations, \
        ContributionMatrixB, EigenvaluesB, OccupationsB = ReadMOComposition(
            Content, LenContent, Type, BasisDimension)

    try:
        from scipy.constants import codata
        HaToeV = codata.value('Hartree energy in eV')
    except ImportError:
        HaToeV = 27.21138602

    Eigenvalues *= HaToeV
    if Type == 'UKS':
        EigenvaluesB *= HaToeV

    fig, ax = plt.subplots()

    PlotParameters = PlotTotalDOS(args, fig, ax, Type, Eigenvalues, Occupations, EigenvaluesB=EigenvaluesB,
                          OccupationsB=OccupationsB)

    SpeciesToPlot = {}

    for key in Species.iterkeys():
        if args.atoms:
            AtomNumber = int(re.split('[a-zA-Z]',key)[0])
            if AtomNumber in args.atoms:
                SpeciesToPlot[key] = Species[key]

    PlotPartialDOS(args, fig, ax, Type, SpeciesToPlot, PlotParameters, Eigenvalues, Occupations, ContributionMatrix,
                   EigenvaluesB=EigenvaluesB, OccupationsB=OccupationsB, ContributionMatrixB=ContributionMatrixB)

    if args.interactive:
        plt.legend()
        plt.show()