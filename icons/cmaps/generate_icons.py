#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Script that generates icons showing the colormaps.
"""

import matplotlib.pyplot as plt
import numpy             as np


if __name__ == '__main__':
    
    cmaps = plt.colormaps()
    
    f      = plt.figure()
    
    for cmap in cmaps:
            
        ax     = f.add_subplot(111)
        f.subplots_adjust(left=0, bottom=0., right=1, top=1)
        
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient, gradient))
        
        ax.imshow(gradient, aspect=10, cmap=cmap)
        ax.set_axis_off()
        
        plt.savefig(f'{cmap}.png', bbox_inches='tight', transparent=True)
        plt.clf()