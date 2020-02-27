import numpy as np
from numpy import load
import scipy.io

data = load('MS1.npz')
lst = data.files

xReceivers = data['anchor'].astype(np.double)
#On ajoute un moins, ils ont inversé les valeurs ces cons...
TDOA = -data['TDOA']
xTag = data['position']

scipy.io.savemat('out.mat', mdict={'xReceivers': xReceivers,'TDOA': TDOA,'xTag': xTag})
