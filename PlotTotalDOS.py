import matplotlib.pyplot as plt
import numpy as np
from functions import Gaussian, dump

def PlotTotalDOS(args, fig, ax, Type, Eigenvalues, Occupations, **kwargs):
    """

    :param fig: Matplotlib figure object
    :param ax: Matplotlib axis object
    :param Type: RKS or UKS
    :param Eigenvalues: if RKS, total, if UKS, spin up MO eigenvalues (in Hartree).
    :param Occupations: if RKS, total, if UKS, spin up MO occupation numbers.
    :param low: low range limit for DOS plot
    :param high: high range limit for DOS plot
    :param smear: gaussian smearing in Hartree
    :param kwargs: if UKS, spin down eigenvalues and occupations go there
    :return: None
    """

    low = args.lowest
    high = args.highest
    smear = args.smear
    PlotParameters = {
        'Domain': None,
        'MoInEnergyRange': None,
        'MoInEnergyRangeB': None
    }

    for key in kwargs:
        if kwargs[key] is not None and key == 'EigenvaluesB':
            EigenvaluesB = kwargs[key]
        elif kwargs[key] is not None and key == 'OccupationsB':
            OccupationsB = kwargs[key]

    MoInEnergyRange = np.where((low < Eigenvalues) & (Eigenvalues < high))[0]
    EnergyRange = Eigenvalues[MoInEnergyRange]
    EnergyLow = EnergyRange[0] - 5 * smear
    EnergyHigh = EnergyRange[-1] + 5 * smear

    if Type == 'UKS':
        MoInEnergyRangeB = np.where((low < EigenvaluesB) & (EigenvaluesB < high))[0]
        EnergyRangeB = EigenvaluesB[MoInEnergyRangeB]
        EnergyLowB = EnergyRangeB[0] - 5 * smear
        EnergyHighB = EnergyRangeB[-1] + 5 * smear
        if EnergyLowB < EnergyLow:
            EnergyLow = EnergyLowB
        if EnergyHighB > EnergyHigh:
            EnergyHigh = EnergyHighB
        PlotParameters['MoInEnergyRangeB'] = np.copy(MoInEnergyRangeB)

    eF = getFermiLevel(Eigenvalues, Occupations) if Type == 'RKS' else getFermiLevel(Eigenvalues, Occupations,
                                                                                     EigenvaluesB=EigenvaluesB,
                                                                                     OccupationsB=OccupationsB)
    Domain = np.linspace(EnergyLow, EnergyHigh, 2001)
    Norm = 2.0 if Type == 'RKS' else 1.0

    PlotParameters['Domain'] = np.copy(Domain)
    PlotParameters['MoInEnergyRange'] = np.copy(MoInEnergyRange)
    TotalDOS = np.zeros_like(Domain)
    for index in MoInEnergyRange:
        TotalDOS += Gaussian(Eigenvalues[index], smear, Norm, Domain)

    if Type == 'UKS':
        TotalDOSB = np.zeros_like(Domain)
        for index in MoInEnergyRangeB:
            TotalDOSB -= Gaussian(EigenvaluesB[index], smear, Norm, Domain)
        plt.plot(Domain, TotalDOSB, '-k', lw=2)
        if args.dump:
            dump(args, '_totalDOS_down.txt',Domain, TotalDOSB)

    plt.plot(Domain, TotalDOS, '-k', lw=2)
    postfix = '_totalDOS.txt' if Type == 'RKS' else '_totalDOS_up.txt'
    dump(args, postfix, Domain, TotalDOS)
    plt.axvline(eF, color='black')
    plt.axhline(0, color='black')

    return PlotParameters

def getFermiLevel(Eigenvalues, Occupations, **kwargs):
    Type = 'RKS'

    for key in kwargs:
        if kwargs[key] is not None and key == 'EigenvaluesB':
            EigenvaluesB = kwargs[key]
        elif kwargs[key] is not None and key == 'OccupationsB':
            OccupationsB = kwargs[key]

    LUMOIndex = None
    for n, (e, f) in enumerate(zip(Eigenvalues, Occupations)):
        if f < 0.5:
            LUMOIndex = n
            HOMOIndex = LUMOIndex - 1
            break

    if Type == 'RKS':
        return 0.5 * (Eigenvalues[HOMOIndex] + Eigenvalues[LUMOIndex])

    elif Type == 'UKS':
        for n, (e, f) in enumerate(zip(EigenvaluesB, OccupationsB)):
            if f < 0.5:
                LUMOBIndex = n
                HOMOBIndex = LUMOBIndex - 1
                break

        HOMO = max(Eigenvalues[HOMOIndex], EigenvaluesB[HOMOBIndex])
        LUMO = min(Eigenvalues[LUMOIndex], EigenvaluesB[LUMOBIndex])
        return 0.5 * (HOMO + LUMO)
