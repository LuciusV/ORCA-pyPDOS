# ORCA-pyPDOS
Partial density-of-states (PDOS) plotter for ORCA (https://orcaforum.cec.mpg.de/) written in python

I. ABOUT
II. REQUIREMENTS
III. FEATURES
IV. FEATURES NOT IMPLEMENTED YET (IN TODO LIST).
V. KNOWN BUGS

I. ABOUT
This python script plots partial density-of-states from ORCA log file.
AUTHOR: Dr. Evgeny V. Tikhonov, Moscow State University, Physics department, email: e.tikhonov@physics.msu.ru

II. REQUIREMENTS
1. ORCA job needs to be completed with keyword
%output print [p_mos] 1 end
2. Python packages required: matplotlib, numpy, scipy, argparse

III. FEATURES
1. Plotting PDOS for each sort of atom
2. Plotting PDOS for every single atom
3. Plotting total DOS

IV. FEATURES NOT IMPLEMENTED YET (IN TODO LIST).
1. Plotting PDOS for list of chosen atoms
2. Interactive and not-interactive mode

V. KNOWN BUGS

(#001) For lines where there is no space between columns, like:
3H 2s -3.975631 0.012837 -2.950634 -0.009080-16.567733 -6.648447
parser fails to split line correctly
WORKAROUND: edit file manually, adding space
FIX IN PROGRESS: fix with regular expression substitution
