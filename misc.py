#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Miscellaneous classes.
"""

import enum
import numpy           as     np
import matplotlib      as     mpl
from   PyQt6.QtWidgets import QToolBar, QComboBox, QWidget

class DummyMouseEvent:
    r'''
    ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
    
    A dummy matplotlib Mouse event that incorporates a subset of its properties useful to mimic mouse events.
    '''
    
    def __init__(self, 
                 canvas : mpl.backend_bases.FigureCanvasBase, 
                 xdata  : int | float, 
                 ydata  : int | float
                ) -> None:
        r'''
        ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
        
        :param canvas: canvas within which the event takes place
        :type canvas: matplotlib.backend_bases.FigureCanvasBase
        :param xdata: x position of the event in figure coordinates
        :type xdata: :python:`int` or :python:`float`
        :param ydata: y position of the event in figure coordinates
        :type ydata: :python:`int` or :python:`float`
        '''
        
        self.canvas = canvas
        self.xdata  = xdata
        self.ydata  = ydata

class Application_states(enum.Enum):
    r'''
    ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
    
    States that the application can be into. This is used to control the interactions with the user (mouse, keyboard, etc.).
    '''
    
    LOCK = enum.auto()

class CustomToolbar(QToolBar):
    r'''
    ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
    
    Custom toolbar with custom widgets.
    '''
    
    def __init__(self, parent: QWidget, root: QWidget, cmap: str, *args, **kwargs) -> None:
        r'''
        ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget holding this widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param root: root widget
        :type root: PyQt6.QtWidgets.QWidget
        :param str cmap: default colormap
        '''
        
        super().__init__('Toolbar', *args, **kwargs)
        
        self.__parent = parent
        self.__root   = root
        
        # Combobox widget containing the list of colormaps
        self.__combobox_cmaps = QComboBox()
        
        # Add matplotlib colormaps to the list of cmaps and set to current cmap
        cmaps = self.root.cmaps_ok
        cmaps.sort()
        
        self.combobox_cmaps.addItems(cmaps)
        self.combobox_cmaps.setCurrentText(cmap)
        
        self.addWidget(self.combobox_cmaps)
        
        ###############################
        #           Signals           #
        ###############################
        
        # When a new cmap is selected, update the image
        self.combobox_cmaps.currentTextChanged.connect(self.update_cmap)
        
        # When the combobox is activated (i.e. clicked with change or not), we give back the focus to the main window
        self.combobox_cmaps.activated.connect(self.root.setFocus)
        
        return
    
    def update_cmap(self, cmap: str) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update cmap and redraw the image when a new cmap is selected in the combobox
        '''
        
        # Stored temporarily the old cmap
        cmap_old = self.root.mpl_im_widget.cmap
        
        # Update cmap and redraw the image
        self.root.mpl_im_widget.cmap = cmap
        self.root.mpl_im_widget.draw()
        
        # Send a status message
        self.root.status_bar.showMessage(f'Colormap was changed from {cmap_old} to {self.root.mpl_im_widget.cmap}.', msecs=3000)
        
        return
    
    ###################################
    #       Getters and setters       #
    ###################################
    
    @property
    def parent(self) -> QWidget:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Parent widget.
        '''
    
        return self.__parent
    
    @property
    def root(self) -> QWidget:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Root widget.
        '''
    
        return self.__root
    
    @property
    def combobox_cmaps(self) -> QComboBox:
        r'''
        ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
        
        Combobox widget containing the list of cmaps from matplotlib.
        '''
        
        return self.__combobox_cmaps

class ArrayList(list):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom list that only allows numpy.ndarray elements of dimensions 2 or 3.
    '''
    
    def __init__(self, *args, **kwargs) -> None:
        
        super().__init__(*args, **kwargs)
        return
    
    def append(self, array: np.ndarray) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Custom append method that checks that the given element is a numpy.ndarray of dimensions 2 or 3.
        
        :param array: 2D or 3D array
        :type array: numpy.ndarray
        
        :raises:
            * :python:`TypeError` if `not isinstance(array, np.ndarray)`
            * :python:`ValueError` if `array.ndim < 2 or array.ndim > 3`
        '''
        
        if not isinstance(array, np.ndarray):
            raise TypeError(f'array to append has type {type(array)} but it must have type numpy.ndarray')

        if array.ndim < 2 or array.ndim > 3:
            raise ValueError(f'array has dimensions {array.ndim} but only data of dimension 2 (images) or 3 (cubes) are allowed.')

        super().append(array)
        return