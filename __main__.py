#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Inspired by PyQubeVis (https://gitlab.lam.fr/bepinat/PyQubeVis/-/tree/master?ref_type=heads) from Epinat Benoit (LAM).
"""

import sys
import signal
import numpy                              as np
import matplotlib.pyplot                  as plt
import matplotlib.figure                  as     mplf
from   PyQt6.QtWidgets                    import QApplication, QMainWindow, QWidget, QVBoxLayout
from   matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class Mpl_im_canvas(FigureCanvas):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom matplotlib canvas that can hold an image.
    
    Heavily inspired by the mplCanvas class in PyQubeVis.
    '''
    
    # Class attribute is the list of allowed cmaps from matplotlib
    __cmaps_ok = plt.colormaps()
    
    def __init__(self, parent: QWidget, cmap: str) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget holding this widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param str cmap: initial colormap used
        
        :raises: 
            * :python:`TypeError` if `not isinstance(cmap, str)`
        '''
        
        #: Parent widget
        self.__parent = parent
        
        self.cmap     = cmap

        # Set the figure and the axis
        self.figure   = mplf.Figure()
        self.__ax     = self.figure.add_subplot(111)
        self.ax.set_axis_off()
        
        # Setup layout of the figure and axes
        self.figure.subplots_adjust(left=0., bottom=0., right=1., top=1.)
        self.ax.axes.tick_params(bottom=False, top=False, left=False, right=False)
        
        # Array containing the data. None means no data has been provided yet
        self.__array  = None
        
        # Artist that contains the image shown
        self.__im_artist = self.ax.imshow(np.full((2, 2), 0), cmap='rainbow')
            
        super().__init__(figure=self.figure)
        
        # Taken from Benoit but commented for now
        '''
        super().setSizePolicy(self, QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
        '''
        super().updateGeometry()
            
        return
    
    ###################################
    #       Getters and setters       #
    ###################################
    
    @property
    def ax(self):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Main axis associated to the figure.
        '''
    
        return self.__ax
    
    @property
    def array(self) -> np.ndarray:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Array that contains the data that are shown as an image.
        '''
        
        return self.__image

    @array.setter
    def array(self, image: np.ndarray | None) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the array that contains the data associated to the image.
        
        :param image: image to show. If None, the image is not initialized.
        :type image: numpy.ndarray
    
        :raises: 
            * :python:`TypeError` if `not isinstance(image, np.ndarray)`
            * :python:`ValueError` if `image.ndim < 2 or image.ndim > 3`
            * :python:`NotImplementedError` if `image.ndim == 3`
        '''
        
        # Checking for the data type of image
        if isinstance(image, None):
            self.__image = None
            
        elif not isinstance(image, np.ndarray):
            raise TypeError(f'Trying to update the image with data of type {type(image)} but only np.ndarray is allowed.')
        
        # Checking for the dimensions of the image
        if image.ndim < 2 or image.ndim > 3:
            raise ValueError(f'Trying to update the image data of shape {image.shape} but data must have 2 or 3 dimensions to be shown.')
            
        elif image.ndim == 3:
            raise NotImplementedError('Showing data cubes is not supported yet...')
        
        self.__image = image
        
        return
    
    @property
    def cmap(self) -> str:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Current colormap used to show the figure.
        '''
        
        return self.__cmap
    
    @cmap.setter
    def cmap(self, cmap: str) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Set the cmap to a new one.
        
        :param str cmap: new colormap for the image
       
        :raises:
            * :python:`TypeError` if `not isinstance(cmap, str)`
            * :python:`ValueError` if the given cmap does not belong to the list of cmaps from matplotlib.pyplot
        '''
        
        if not isinstance(cmap, str):
            raise TypeError(f'cmap is of type {type(cmap)} but it must be of type str.')
        
        if cmap not in Mpl_im_canvas.__cmaps_ok:
            raise ValueError(f'cmap is {cmap} which does not belong to the following list of cmaps from matplotlib: {Mpl_im_canvas.__cmaps_ok}.')
        
        self.__cmap   = cmap

        return
    
    ##################################################
    #        Updating properties of the image        #
    ##################################################
    
    def update_image(self, image: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param image: image to show
        :type image: numpy.ndarray
        '''
        
        # Update the array
        self.array = image
        
        # If first instantiation of the artist
        if self.__im_artist is None:
            self.__im_artist = self.ax.imshow(self.array, cmap=self.cmap)
            
        # Otherwise, just update the data of the artist
        else:
            self.im_artist.set_data(self.array)
          
        # Apply changes
        self.updateGeometry()
        
        return
            
class Mpl_image_widget(QWidget):
    
    def __init__(self, parent):
        
        self.parent = parent
        
        super().__init__(self, self.parent)
        
        # Setup custom canvas
        self.canvas = Mpl_im_canvas()
        
        self.vbl = QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
        

class Window(QMainWindow):
    
    def __init__(self, *arg, **kwargs):
        
        super().__init__(*arg, **kwargs)
        
        # Main frame that serves as central widget
        #mpl_im_widget = Mpl_image_widget(self)
        mpl_im_widget = Mpl_im_canvas(self, 'rainbow')
        self.setCentralWidget(mpl_im_widget)
        
        # Add
        
        return

def main(argv):
    
    # Main application and window
    app = QApplication(sys.argv)
    
    # Handling SIGINT signal from terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    main_window = Window()
    main_window.show()
    
    # Add the possibility to close the window through the terminal
    sys.exit(app.exec())

if __name__ == '__main__':
    main(sys.argv[1:])