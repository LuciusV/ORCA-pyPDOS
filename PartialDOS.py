from ReadOrcaLog import ReadOrcaLog
from ReadMOComposition import ReadMOComposition

if __name__ == '__main__':

    OrcaLog = '/var/tmp/O2.out'
    f = open(OrcaLog, 'rb')
    Type, NumberOfElectrons, BasisDimension, Content, LenContent = ReadOrcaLog(f)
    f.close()
    Target = 'MOLECULAR ORBITALS'
    if Type == 'RKS':
        Species, ContributionMatrix, Eigenvalues, Occupations = ReadMOComposition(Content, LenContent, Type, BasisDimension)
    elif Type == 'UKS':
        Species, ContributionMatrixA, EigenvaluesA, OccupationsA, \
        ContributionMatrixB, EigenvaluesB, OccupationsB = ReadMOComposition(
            Content, LenContent, Type, BasisDimension)

    print('')