import numpy as np

def Gaussian(Center, Smear, Norm, Domain):
    return Norm * np.exp( - 0.5 * (Domain - Center) ** 2 / Smear ** 2 )

def dump(args, suffix, x, y):
    fname = args.filename.split('.out')[0] + suffix
    np.savetxt(fname, np.transpose([x,y]))