# This test should be placed in some 'test' folder, but currently lies in the root folder in order to easily import gate.protonct and pctpairprotons.
# It should be moved once the testing mechanics are finalized.

import os
import itk
import numpy as np

from gate.protonct import protonct
from pctpairprotons import pctpairprotons

test_folder = '/tmp/pctPairProtons'

protonct(
    output=test_folder,
    projections=1,
    protons_per_projection=100,
    seed=123
)

pctpairprotons(
    input_in=os.path.join(test_folder, 'PhaseSpaceIn.root'),
    input_out=os.path.join(test_folder, 'PhaseSpaceOut.root'),
    output=os.path.join(test_folder, 'pairs.mhd'),
    plane_in=-110,
    plane_out=110,
    psin='PhaseSpaceIn',
    psout='PhaseSpaceOut'
)

im = itk.imread(os.path.join(test_folder, 'pairs0000.mhd'))
assert np.asarray(im).shape == (97, 5, 3)
