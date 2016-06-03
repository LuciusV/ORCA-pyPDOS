def ReadOrcaLog(fileobject):
    """
    :param fileobject: open file object of Orca job log
    :return Type: RKS or UKS
    :return int(NumberOfElectrons): Number of electrons
    :return int(BasisDimension): Basis dimension (equals to number of molecular orbitals).
    :return Content: list of lines in file
    :return LastLine: number of lines in file
    """
    Content = fileobject.readlines()
    Type = None
    NumberOfElectrons = None
    BasisDimension = None
    index = 0
    LastLine = len(Content)
    while not (Type and NumberOfElectrons and BasisDimension):
        if 'Number of Electrons' in Content[index]:
            NumberOfElectrons = Content[index].split()[-1]
        elif 'Basis Dimension' in Content[index]:
            BasisDimension = Content[index].split()[-1]
        elif 'Hartree-Fock type' in Content[index]:
            if 'UHF' in Content[index]:
                Type = 'UKS'
            elif 'RHF' in Content[index]:
                Type = 'RKS'
        index += 1
        if index == LastLine:
            raise TypeError('Orca log file incomplete or wrong!')

    return Type, int(NumberOfElectrons), int(BasisDimension), Content, LastLine