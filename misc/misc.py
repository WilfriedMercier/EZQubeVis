#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Miscellaneous classes.
"""

import enum
import numpy                              as     np
import matplotlib                         as     mpl
import os.path                            as     opath
from   functools                          import partialmethod
from   PyQt6.QtWidgets                    import QToolBar, QComboBox, QWidget, QLabel
from   PyQt6.QtGui                        import QIcon
from   PyQt6.QtCore                       import QSize, Qt, QEvent

class Base_widget_skeleton:
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    A simple mother class common to all widgets used in this application that stores the parent and root widgets as private variables and define getters only.
    '''
    
    def __init__(self, parent: QWidget, root: QWidget) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param root: root widget
        :type root: PyQt6.QtWidgets.QWidget
        
        :raises TypeError: if one of the following
            * :python:`not isinstance(parent, PyQt6.QtWidgets.QWidget)`
            * :python:`not isinstance(root, PyQt6.QtWidgets.QWidget)`
        '''
        
        if not isinstance(parent, QWidget):
            raise TypeError(f'Parent widget has type {type(parent)} but it should be of type PyQt6.QtWidgets.QWidget.')
                
        if not isinstance(root, QWidget):
            raise TypeError(f'Root widget has type {type(root)} but it should be of type PyQt6.QtWidgets.QWidget.')
        
        self.__parent = parent
        self.__root   = root
        
        return
        
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
    
class Dummy_mouse_event:
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
    
    # Lock state: the highlight rectangle is locked on the image and the spectrum is not updated when the mouse moves
    LOCK = enum.auto()
    
    # Mask state: locking is disabled, when the mouse is pressed and held the user can draw a rectangle to mask pixels
    MASK = enum.auto()
    
    # Mask activated: only when the mask is activated, does the masking work. The activation happens when the mouse is kept pressed
    MASK_ON = enum.auto()
    
class Data_types(enum.Enum):
    r'''
    ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
    
    Different types of input data that the program can handle.
    '''
    
    # Flag used for a 2D array corresponding to an image
    IMAGE = enum.auto()
    
    # Flag used for a 3D array corresponding to a cube
    CUBE  = enum.auto()

class Combobox_cmaps(Base_widget_skeleton, QComboBox):
    r'''
    ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
    
    A custom combobox made to show colormaps.
    '''
    
    def __init__(self, parent: QWidget, root: QWidget, cmap: str, *args, **kwargs):
        r'''
        ..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>
        
        :param str cmap: default colormap
        '''
        
        Base_widget_skeleton.__init__(self, parent, root)
        QComboBox.__init__(self, *args, **kwargs)
        
        # Cmap selected kept in memory to allow cmap preview on the fly
        self.__cmap_selected = cmap
        
        # Add matplotlib colormaps to the list of cmaps and set to current cmap
        self.addItems(self.root.cmaps_ok)
        self.setCurrentText(cmap)
        self.setFrame(False)
        
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        
        ###############################
        #           Signals           #
        ###############################
        
        # When a new cmap is selected, update the image
        self.currentTextChanged.connect(self.update_cmap)
        
        self.activated.connect(self.root.setFocus)
        
        # When a cmap is highlighted, we show how it looks. If not selected, we go back to the original cmap
        self.highlighted.connect(lambda index: self.preview_cmap(self.root.cmaps_ok[index]))
        
        # Event filter is mandatory to handle a custom escape key event
        self.view().installEventFilter(self)
    
    #################################
    #       Signals and slots       #
    #################################
    
    def eventFilter(self, obj, event) -> bool:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Custom event filter that handles the escape key in the combobox.
        '''
        
        if event.type() == QEvent.Type.KeyPress and event.key() == Qt.Key.Key_Escape:
            
            # Reinitialize the cmap first
            self.update_cmap(self.__cmap_selected)
            
            # Do what the escape key is supposed to do
            super().eventFilter(obj, event)
        
        return super().eventFilter(obj, event)
        
    def update_cmap(self, cmap: str, preview: bool = False) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update cmap and redraw the image when a new cmap is selected in the combobox
        
        .. note::
            
            A preview will not store the cmap name. Thus, when the escape key is pressed, the cmap will be restored to the value stored in **self.__cmap_selected**
        
        :param str cmap: new colormap
         
        Keyword arguments
        -----------------
        
        :param bool preview: whether this is a preview or a definite update of the cmap
        '''
        
        # Store the new cmap
        if not preview:
            self.__cmap_selected = cmap
        
        # Store temporarily the old cmap
        cmap_old = self.root.mpl_im_widget.cmap
        
        # Update cmap and redraw the image
        self.root.mpl_im_widget.cmap = cmap
        self.root.mpl_im_widget.draw()
        
        # Send a status message
        self.root.status_bar.showMessage(f'Colormap was changed from {cmap_old} to {self.root.mpl_im_widget.cmap}.', msecs=3000)
        
        return
    
    # Method that updates the colormap of the image but in preview mode
    preview_cmap = partialmethod(update_cmap, preview=True)

class Custom_toolbar(Base_widget_skeleton, QToolBar):
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
        
        Base_widget_skeleton.__init__(self, parent, root)
        QToolBar.__init__(self, 'Toolbar', *args, **kwargs)
        
        
        # Label before Combobox
        self.__label_cmaps   = QLabel('Colormap ')
        
        self.addWidget(self.__label_cmaps)
        
        # Combobox widget containing the list of colormaps
        self.__combobox_cmaps = Combobox_cmaps(self, self.root, cmap)
        
        for pos, cmap in enumerate(self.root.cmaps_ok):
            self.combobox_cmaps.setItemIcon(pos, QIcon(opath.join(f'icons/cmaps/{cmap}.png')))
        
        self.combobox_cmaps.setIconSize(QSize(100, 10))
        
        self.addWidget(self.combobox_cmaps)
        
        return
    
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