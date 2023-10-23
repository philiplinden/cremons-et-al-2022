from matplotlib import pyplot as plt
import numpy as np
from pprint import pprint as print

from simulator import Hapke, get_sample_grid


if __name__ == '__main__':
    WLS = get_sample_grid()
    Refl = np.array([np.random.randn()*x for x in WLS])
    R = Hapke(Refl, WLS)
    print(WLS)
    print(Refl)
    print(R)
    plt.figure()
    plt.plot(WLS,Refl,WLS,R)
    plt.show()
