from functions import Gaussian, dump
import numpy as np
import matplotlib.pyplot as plt

def PlotPartialDOS(args, fig, ax, Type, SpeciesToPlot, PlotParameters, Eigenvalues, Occupations, ContributionMatrix, **kwargs):

    Domain = PlotParameters['Domain']
    MoInEnergyRange = PlotParameters['MoInEnergyRange']
    MoInEnergyRangeB = PlotParameters['MoInEnergyRangeB']

    for key in kwargs:
        if kwargs[key] is not None and key == 'EigenvaluesB':
            EigenvaluesB = kwargs[key]
        elif kwargs[key] is not None and key == 'OccupationsB':
            OccupationsB = kwargs[key]
        elif kwargs[key] is not None and key == 'ContributionMatrixB':
            ContributionMatrixB = kwargs[key]

    SumOfselected = np.zeros_like(Domain)
    if Type == 'UKS':
        SumOfselectedB = np.zeros_like(Domain)

    for key in SpeciesToPlot:
        AtomDos = np.zeros_like(Domain)
        for index in MoInEnergyRange:
            norm = np.sum(ContributionMatrix[index,:]**2)
            if Type == 'RKS':
                norm *= 0.5
            weight = np.sum(ContributionMatrix[index,range(SpeciesToPlot[key][0], SpeciesToPlot[key][1]+1)]**2) / norm
            AtomDos += Gaussian(Eigenvalues[index],args.smear, weight, Domain)
        if Type == 'UKS':
            AtomDosB = np.zeros_like(Domain)
            for index in MoInEnergyRangeB:
                norm = np.sum(ContributionMatrixB[index, :] ** 2)
                weight = np.sum(ContributionMatrixB[index, range(SpeciesToPlot[key][0], SpeciesToPlot[key][1] + 1)] ** 2) / norm
                AtomDosB -= Gaussian(EigenvaluesB[index], args.smear, weight, Domain)
        postfix = '' if Type == 'RKS' else ' up'
        SumOfselected += AtomDos
        if args.unique:
            dump(args, '_pdos_' + key + postfix.replace(' ','_'), Domain, AtomDos)
            plt.plot(Domain, AtomDos, label=key+postfix)
        if Type == 'UKS':
            SumOfselectedB += AtomDosB
            if args.unique:
                dump(args,'_pdos_'+key+'_down.txt',Domain, AtomDosB)
                plt.plot(Domain, AtomDosB, label=key+' down')

    plt.plot(Domain, SumOfselected, lw=2, label='full of sel. atoms' + postfix)
    if Type == 'UKS':
        postfix = postfix.replace(' ','_')
    dump(args, '_pdos'+postfix+'.txt',Domain, SumOfselected)
    if Type == 'UKS':
        plt.plot(Domain, SumOfselectedB, lw=2, label='full of sel. atoms down')
        if args.dump:
            dump(args, '_pdos_down.txt',Domain,SumOfselectedB)