import numpy as np

def ReadMOComposition(Content, LenContent, Type, MONumber):
    """
    This function makes arrays of Eigenvalues, Occupations and matrix of basis set orbitals contribution into MOs.
    [i,j] is contribution into i-th molecular orbital from j-th basis set function.
    :param Content: ORCA Log file as list of strings
    :param LenContent: Length of that list
    :param Type: RKS or UKS
    :param MONumber: Basis set dimension (equal to number of MOs)
    :return Species: list of different atoms
    :return ContributionMatrix: if RKS, total, if UKS, spin up contribution matrix.
    :return Eigenvalues: if RKS, total, if UKS, spin up MO eigenvalues (in Hartree).
    :return Occupations: if RKS, total, if UKS, spin up MO occupation numbers.
    :return ContributionMatrixB: Only if UKS, spin down contribution matrix.
    :return EigenvaluesB: Only if UKS, spin down MO eigenvalues (in Hartree).
    :return OccupationsB: Only if UKS, spin down MO occupation numbers.
    """

    offset = None
    for n, line in enumerate(reversed(Content)):
        if 'MOLECULAR ORBITALS' in line:
            offset = LenContent - n + 1
            break

    if not offset:
        raise TypeError(
            'Molecular orbitals not found. Please rerun ORCA job with option "%output print [p_mos] 1 end".')

    Species = []
    ContributionMatrix = np.zeros((MONumber, MONumber))
    Eigenvalues = np.zeros(MONumber)
    Occupations = np.zeros(MONumber)

    '''
    Here we read MO coefficients from columns of 6.
    Last page of columns can be less than 6.
    Page is from 6 MO numbers until last coefficients before next page,
    so its total length is basis set dimension (variable MONumber) + 4.
    '''

    EndPage = MONumber // 6

    Columns = 6

    if MONumber % 6 != 0:
        EndPage += 1

    for index in xrange(EndPage * (MONumber + 4)):
        Line = Content[offset + index][:-1]
        LineOnPage = index % (MONumber + 4)
        Page = index // (MONumber + 4)
        if Page + 1 == EndPage:
            Columns = MONumber % 6
        if LineOnPage == 1:  # Eigenvalues
            # print Page, EndPage, [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in range(Columns)]
            # print [float(x) for x in [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in range(Columns)]]
            Eigenvalues[Page * 6:Page * 6 + Columns] = [float(E) for E in
                                                        [Line[shift * 10 + 17: shift * 10 + 17 + 9] for shift in
                                                         range(Columns)]]
        elif LineOnPage == 2:  # Occupations
            # print [Line[shift * 10 + 19: shift * 10 + 19 + 7] for shift in range(Columns)]
            # print [float(x) for x in [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in range(Columns)]]
            Occupations[Page * 6:Page * 6 + Columns] = [float(E) for E in
                                                        [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in
                                                         range(Columns)]]
        elif LineOnPage >= 4:
            ContributionMatrix[Page * 6:Page * 6 + Columns, (LineOnPage - 4 + Page * MONumber) % MONumber] = \
                [float(E) for E in [Line[shift * 10 + 17: shift * 10 + 17 + 9] for shift in range(Columns)]]
            if Page == 0:
                Atom = str(Line[:5]).strip(' ')
                if Atom not in Species:
                    Species.append(Atom)

    if type == 'RKS':
        return Species, ContributionMatrix, Eigenvalues, Occupations

    if Type == 'UKS':
        # repeat again for spin down! Sorry for copy-paste...
        ContributionMatrixB = np.zeros((MONumber, MONumber))
        EigenvaluesB = np.zeros(MONumber)
        OccupationsB = np.zeros(MONumber)

        offset += (MONumber + 4) * EndPage + 1
        Columns = 6

        for index in xrange(EndPage * (MONumber + 4)):
            Line = Content[offset + index][:-1]
            LineOnPage = index % (MONumber + 4)
            Page = index // (MONumber + 4)
            if Page + 1 == EndPage:
                Columns = MONumber % 6
            if LineOnPage == 1:  # Eigenvalues
                # print Page, EndPage, [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in range(Columns)]
                # print [float(x) for x in [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in range(Columns)]]
                EigenvaluesB[Page * 6:Page * 6 + Columns] = [float(E) for E in
                                                             [Line[shift * 10 + 17: shift * 10 + 17 + 9] for shift in
                                                              range(Columns)]]
            elif LineOnPage == 2:  # Occupations
                # print [Line[shift * 10 + 19: shift * 10 + 19 + 7] for shift in range(Columns)]
                # print [float(x) for x in [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in range(Columns)]]
                OccupationsB[Page * 6:Page * 6 + Columns] = [float(E) for E in
                                                             [Line[shift * 10 + 18: shift * 10 + 18 + 8] for shift in
                                                              range(Columns)]]
            elif LineOnPage >= 4:
                ContributionMatrixB[Page * 6:Page * 6 + Columns, (LineOnPage - 4 + Page * MONumber) % MONumber] = \
                    [float(E) for E in [Line[shift * 10 + 17: shift * 10 + 17 + 9] for shift in range(Columns)]]

        return Species, ContributionMatrix, Eigenvalues, Occupations, ContributionMatrixB, EigenvaluesB, OccupationsB

