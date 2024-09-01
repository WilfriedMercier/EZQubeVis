#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Custom classes built around the FigureCanvas class from matplotlib.
"""

import matplotlib
import numpy                              as     np
import matplotlib.figure                  as     mplf
from   PyQt6.QtWidgets                    import QWidget, QDockWidget
from   PyQt6.QtCore                       import Qt
from   matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Custom import
from   .misc                              import (Application_states, Dummy_mouse_event, 
                                                  Base_widget_skeleton
                                                 )

class Mpl_spectrum_canvas(FigureCanvas):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom matplotlib canvas that can hold a spectrum.
    
    Heavily inspired by the mplCanvas class in PyQubeVis.
    '''
    
    def __init__(self,
                 parent           : QWidget,
                 root             : QWidget,
                 wavelength       : np.ndarray | None = None,
                 wavelength_model : np.ndarray | None = None,
                 wavelength_unit  : str        | None = None,
                 spec_pos         : int        | None = None
                ) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget holding this widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param root: root widget
        :type root: PyQt6.QtWidgets.QWidget
        
        Keyword arguments
        -----------------
        
        :param wavelength: wavelength for the data spectrum
        :type wavelength: numpy.ndarray or :python:`None`
        :param wavelength_model: wavelength for the model spectrum
        :type wavelength_model: numpy.ndarray or :python:`None`
        :param wavelength_unit: unit of wavelength. If None, no wavelength is shown.
        :type wavelength_unit: :python:`str` or :python:`None`
        :param spec_pos: initial position in the spectrum that is shown on the image if a cube is shown. If a 2D image is shown instead, this parameter is None.
        :type spec_pos: :python:`int` or :python:`None`
        
        :raise: TypeError if
            * :python:`not isinstance(parent, QWidget)`
            * :python:`wavelength is not None and not isinstance(wavelength, np.ndarray)`
            * :python:`wavelength_model is not None and not isinstance(wavelength_model, np.ndarray)`
            * :python:`wavelength_unit is not None and not isinstance(wavelength_unit, str)`
            * 
        '''
        
        if not isinstance(parent, QWidget):
            raise TypeError(f'Parent widget has type {type(parent)} but it should be of type QWidget.')
            
        if not isinstance(root, QWidget):
            raise TypeError(f'Root widget has type {type(root)} but it should be of type QWidget.')
            
        if wavelength is not None and not isinstance(wavelength, np.ndarray):
            raise TypeError(f'wavelength range has type {type(wavelength)} but it should be None or numpy.ndarray.')
            
        if wavelength_model is not None and not isinstance(wavelength_model, np.ndarray):
            raise TypeError(f'model wavelength range has type {type(wavelength_model)} but it should be None or numpy.ndarray.')
        
        if wavelength_unit is not None and not isinstance(wavelength_unit, str):
            raise TypeError(f'wavelength unit has type {type(wavelength_unit)} but it should be None or str.')
        
        if spec_pos is not None and not isinstance(spec_pos, int):
            raise TypeError(f'spectrum position has type {type(spec_pos)} but it should be None or int.')
        
        # Parent widget
        self.__parent = parent

        # Root widget
        self.__root   = root
        
        # Artist that holds the spectrum
        # spec_artist is set private and has neither getter nor setter
        # Initialization to None before creating a first artist in self.update_spectrum when a spectrum is provided
        self.__spec_artist = None
        
        # Artist that holds the model spectrum
        # spec_model_artist is set private and has neither getter nor setter
        # Initialization to None before creating a first artist in self.update_spectrum when a spectrum is provided
        self.__spec_model_artist = None
        
        # Save the wavelength ranges if provided
        self.__wavelength       = wavelength
        self.__wavelength_model = wavelength_model
        
        # Artist that adds a vertical line to show the position in the cube if no input image is provided
        self.__spec_position_line  = None
        
        ########################################
        #                Figure                #
        ########################################
        
        # Note that the figure cannot be private because of the parent class that requires it to be public
        self.figure   = mplf.Figure()
        
        super().__init__(figure=self.figure)
        
        # Axis is set private and has no setter
        self.__ax     = self.figure.add_subplot(111)
        
        # Setup layout of the figure and axes
        self.figure.subplots_adjust(left=0.01, right=0.99, top=1., bottom=0.15)
        self.ax.axes.tick_params(direction='in', bottom=False, top=False, left=False, right=False, labelleft=False, labelright=False, labeltop=False, labelbottom=False)
        self.ax.spines[['left', 'right', 'top', 'bottom']].set_visible(False)
        
        return
    
    ###################################
    #       Getters and setters       #
    ###################################
    
    @property
    def parent(self):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Parent widget.
        '''
    
        return self.__parent
    
    @property
    def root(self):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Root widget.
        '''
    
        return self.__root
    
    @property
    def ax(self):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Main axis associated to the figure.
        '''
    
        return self.__ax
    
    @property
    def spectrum(self) -> np.ndarray:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Array that contains the y-data of the spectrum that is shown in the plot.
        '''
        
        return self.__spectrum

    @spectrum.setter
    def spectrum(self, spectrum: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the array that contains the y-data of the spectrum associated to the plot.
        
        :param spectrum: spectrum to show
        :type spectrum: numpy.ndarray
    
        :raises: 
            * :python:`TypeError` if `not isinstance(spectrum, np.ndarray)`
            * :python:`ValueError` if `spectrum.ndim != 1`
        '''
            
        # Otherwise we check that this is a numpy.ndarray first
        if not isinstance(spectrum, np.ndarray):
            raise TypeError(f'Trying to update the spectrum with data of type {type(spectrum)} but only np.ndarray is allowed.')
        
        # Checking for the dimensions of the image
        if spectrum.ndim != 1:
            raise ValueError(f'Trying to update the spectrum with data of shape {spectrum.shape} but data must have 1 dimension only.')

        # If everything is fine, we just store the data
        self.__spectrum = spectrum
        
        return
    
    @property
    def model_spectrum(self) -> np.ndarray:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Array that contains the y-data of the model that is shown in the plot.
        '''
        
        return self.__model_spectrum

    @model_spectrum.setter
    def model_spectrum(self, spectrum: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the array that contains the y-data of the model spectrum associated to the plot.
        
        :param spectrum: spectrum to show
        :type spectrum: numpy.ndarray
    
        :raises: 
            * :python:`TypeError` if `not isinstance(spectrum, np.ndarray)`
            * :python:`ValueError` if `spectrum.ndim != 1`
        '''
            
        # Otherwise we check that this is a numpy.ndarray first
        if not isinstance(spectrum, np.ndarray):
            raise TypeError(f'Trying to update the spectrum with data of type {type(spectrum)} but only np.ndarray is allowed.')
        
        # Checking for the dimensions of the image
        if spectrum.ndim != 1:
            raise ValueError(f'Trying to update the spectrum with data of shape {spectrum.shape} but data must have 1 dimension only.')

        # If everything is fine, we just store the data
        self.__model_spectrum = spectrum
        
        return
    
    @property
    def wavelength(self) -> np.ndarray:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Array that contains wavelengths corresponding to the data points.
        '''
        
        return self.__wavelength
    
    @property
    def wavelength_model(self) -> np.ndarray:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Array that contains wavelengths corresponding to the model points.
        '''
        
        return self.__wavelength_model
    
    @property
    def cube_position_line(self) -> matplotlib.lines.Line2D | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Line that shows the position in the cube.
        '''
    
        return self.__cube_position_line
    
    #####################################################
    #        Updating properties of the spectrum        #
    #####################################################
    
    def update_spectrum(self, spectrum: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the figure with the given spectrum.
        
        :param spectrum: 1D spectrum to show
        :type spectrum: numpy.ndarray
        '''
        
        # Update the y-data of the spectrum
        self.spectrum = spectrum
        
        # If first instantiation of the artist
        if self.__spec_artist is None:
            
            # If a wavelength range is provided
            if self.wavelength is not None:
                
                self.__spec_artist = self.ax.plot(self.wavelength, spectrum, color='k', lw=1.5)
                
                # Update the layout of the figure
                self.ax.spines['bottom'].set_visible(True)
                self.ax.tick_params(bottom=True, labelbottom=True)
                self.ax.set_xlabel('Wavelength')
                
            # If no wavelength range is provided
            else:
                self.__spec_artist = self.ax.plot(spectrum, color='k', lw=1.5)
                
            self.__spec_artist[0].set_drawstyle('steps')
                
            
        # Otherwise, just update the y-data of the artist
        else:
            
            self.__spec_artist[0].set_ydata(spectrum)
            self.ax.set_ylim([0.9*np.nanmin(spectrum), 1.1*np.nanmax(spectrum)])
        
        return
    
    def update_model_spectrum(self, spectrum: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the figure with the given model spectrum.
        
        :param spectrum: 1D model spectrum to show
        :type spectrum: numpy.ndarray
        '''
        
        # Update the y-data of the spectrum
        self.model_spectrum = spectrum
        
        # If first instantiation of the artist
        if self.__spec_model_artist is None:
            
            if self.wavelength_model is not None:
                
                self.__spec_model_artist = self.ax.plot(self.wavelength_model, self.spectrum, color='firebrick', lw=1)
                
                # Update the layout of the figure
                self.ax.spines['bottom'].set_visible(True)
                self.ax.tick_params(bottom=True, labelbottom=True)
                self.ax.set_xlabel('Wavelength')
                
            else:
                self.__spec_model_artist = self.ax.plot(self.spectrum, color='firebrick', lw=1)
                
            self.__spec_model_artist[0].set_drawstyle('steps')
            
        # Otherwise, just update the y-data of the artist
        else:
            self.__spec_model_artist[0].set_ydata(self.model_spectrum)
        
        return
    
    def update_spec_line_position(self) -> None:
        
        # Create line on firt call
        if self.__spec_position_line is None:
            self.__spec_position_line = self.ax.axvline(self.wavelength[self.root.cube_pos], color='k', ls='--')
        else:
            self.__spec_position_line.set_xdata([self.wavelength[self.root.cube_pos]]*2)
        
        return
    
class Dock_widget_spectrum(QDockWidget):
    
    def __init__(self, 
                 parent           : QWidget, 
                 root             : QWidget,
                 wavelength       : np.ndarray | None = None,
                 wavelength_model : np.ndarray | None = None,
                 spec_pos         : int        | None = None
                ):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param root: root widget
        :type root: PyQt6.QtWidgets.QWidget
        
        Keyword arguments
        -----------------
        
        :param wavelength: wavelength for the data spectrum
        :type wavelength: numpy.ndarray or :python:`None`
        :param wavelength_model: wavelength for the model spectrum
        :type wavelength_model: numpy.ndarray or :python:`None`
        :param spec_pos: initial position in the spectrum that is shown on the image if a cube is shown. If a 2D image is shown instead, this parameter is None.
        :type spec_pos: :python:`int` or :python:`None`
        '''
        
        self.__parent      = parent
        self.__root        = root
        
        super().__init__()
        
        # Canvas holding the spectrum
        self.__spec_canvas = Mpl_spectrum_canvas(self, self.root,
                                                 wavelength       = wavelength, 
                                                 wavelength_model = wavelength_model,
                                                 spec_pos         = spec_pos
                                                )
        
        self.setWidget(self.spec_canvas)
        
        # Disable the possibility to close the dock
        self.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetMovable)
        
    @property
    def parent(self) -> QWidget:
        
        return self.__parent
    
    @property
    def root(self) -> QWidget:
        
        return self.__root
        
    @property
    def spec_canvas(self) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Canvas holding the spectrum
        '''
        
        return self.__spec_canvas