#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Custom classes built around the FigureCanvas class from matplotlib.
"""

import time
import matplotlib
import numpy                              as     np
import matplotlib.pyplot                  as     plt
import matplotlib.figure                  as     mplf
from   PyQt6.QtWidgets                    import QWidget, QDockWidget
from   matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class Mpl_im_canvas(FigureCanvas):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom matplotlib canvas that can hold an image.
    
    Heavily inspired by the mplCanvas class in PyQubeVis.
    '''
    
    # Class attribute is the list of allowed cmaps from matplotlib
    __cmaps_ok = plt.colormaps()
    
    def __init__(self, 
                 parent : QWidget, 
                 cmap   : str,
                ) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget holding this widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param str cmap: initial colormap used
        '''
        
        #: Parent widget
        self.__parent = parent
        
        # Using the setter defined below to automatically perform the checks
        self.cmap     = cmap

        # Set the figure and the axis
        # Note that figure cannot be private because of the parent class that requires it to be public
        self.figure   = mplf.Figure()
        
        super().__init__(figure=self.figure)
        
        # Axis is set private and has no setter
        self.__ax     = self.figure.add_subplot(111)
        self.ax.set_axis_off()
        
        # Setup layout of the figure and axes
        self.figure.subplots_adjust(left=0., bottom=0., right=1., top=1.)
        self.ax.axes.tick_params(bottom=False, top=False, left=False, right=False)
        
        # im_artist is set private and has neither getter nor setter
        # Initialization to None before creating a first artist in self.update_image when an image is provided
        self.__im_artist = None
        
        ######################################
        #           Show the image           #
        ######################################
        
        if self.__parent.image is not None:
            image = self.__parent.image
        elif self.__parent.cube is not None:
            image = self.__parent.cube[0, :, :]
        elif self.__parent.cube_model is not None:
            image = self.__parent.cube_model[0, :, :]
        else:
            raise ValueError('No file provided, cannot show an image.')
            
        self.update_image(image)
        self.draw()
        
        # Update the vmin and vmax values of the image
        if self.array is None:
            self.vmin = None
            self.vmax = None
        else:
            self.vmin = np.nanquantile(self.array, 0.1)
            self.vmax = np.nanquantile(self.array, 0.9)
        
        # Plot artist that adds a square around the currently selected pixel
        self.__highlight_rect = None
        
        #######################################
        #           Events handling           #
        #######################################
        
        self.mpl_connect('motion_notify_event', self.mouse_move)
            
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
    def highlight_rect(self) -> matplotlib.lines.Line2D | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Line2D object used to highlight the currently selected pixel.
        '''
    
        return self.__highlight_rect
    
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
        
        return self.__array

    @array.setter
    def array(self, image: np.ndarray | None) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the array that contains the data associated to the image.
        
        :param image: image to show. If None, the image is not initialized.
        :type image: numpy.ndarray
    
        :raises: 
            * :python:`TypeError` if `not isinstance(image, np.ndarray)`
            * :python:`ValueError` if `image.ndim != 2`
        '''
            
        # Otherwise we check that this is a numpy.ndarray first
        if not isinstance(image, np.ndarray):
            raise TypeError(f'Trying to update the image with data of type {type(image)} but only np.ndarray is allowed.')
        
        # Checking for the dimensions of the image
        if image.ndim != 2:
            raise ValueError(f'Trying to update the image data of shape {image.shape} but data must be 2D.')
       
        # If everything is fine, we just store the data
        self.__array = image
        
        # Update vmin and vmax
        self.vmin    = np.nanquantile(self.array, 0.2)
        self.vmax    = np.nanquantile(self.array, 0.8)
        
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
            self.__im_artist = self.ax.imshow(self.array, 
                                              cmap   = self.cmap, 
                                              vmin   = self.vmin, 
                                              vmax   = self.vmax,
                                              origin = 'lower'
                                             )
            
        # Otherwise, just update the data of the artist
        else:
            self.__im_artist.set_data(self.array)
        
        return
    
    #########################################
    #           Mouse interaction           #
    #########################################
    
    def mouse_move(self, event) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Function called when the mouse is moved. It updates the spectrum at the current cursor location.
        '''
        
        # If None, we are outside and we do not update
        if event.xdata is None or event.ydata is None:
            return
        
        xpos = int(np.round(event.xdata))
        ypos = int(np.round(event.ydata))
        
        ####################################
        #           Image update           #
        ####################################
        
        xvals = [xpos - 0.5, xpos + 0.5, xpos + 0.5, xpos - 0.5, xpos - 0.5]
        yvals = [ypos - 0.5, ypos - 0.5, ypos + 0.5, ypos + 0.5, ypos - 0.5]
        
        if self.highlight_rect is None:
            self.__highlight_rect, = self.ax.plot(xvals, yvals, lw=3, color='k')
        else:
            self.highlight_rect.set_xdata(xvals)
            self.highlight_rect.set_ydata(yvals)
        
        #######################################
        #           Spectrum update           #
        #######################################
        
        # Canvas holding the spectrum to be updated
        spec_canvas = self.parent.bottom_dock.spec_canvas
        
        # Only do the spectrum update if a data cube is provided
        if self.parent.cube is not None:
            
            # Update the data spectrum
            spec_canvas.update_spectrum(self.parent.cube[:, ypos, xpos])
            
        # Only do the spectrum update if a model cube is provided
        if self.parent.cube_model is not None:
            spec_canvas.update_model_spectrum(self.parent.cube_model[:, ypos, xpos])
        
        # Apply changes to the spectrum
        spec_canvas.draw()
        
        # Apply changes to the image
        event.canvas.draw()
        
        return
    
class Mpl_spectrum_canvas(FigureCanvas):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom matplotlib canvas that can hold a spectrum.
    
    Heavily inspired by the mplCanvas class in PyQubeVis.
    '''
    
    def __init__(self,
                 parent   : QWidget,
                 wavelength       : np.ndarray | None = None,
                 wavelength_model : np.ndarray | None = None,
                 wavelength_unit  : str | None        = None
                ) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget holding this widget
        :type parent: PyQt6.QtWidgets.QWidget
        
        Keyword arguments
        -----------------
        
        :param wavelength: wavelength for the data spectrum
        :type wavelength: numpy.ndarray or :python:`None`
        :param wavelength_model: wavelength for the model spectrum
        :type wavelength_model: numpy.ndarray or :python:`None`
        :param wavelength_unit: unit of wavelength. If None, no wavelength is shown.
        :type wavelength_unit: :python:`str` or :python:`None`
        
        :raise: TypeError if
            * :python:`not isinstance(parent, QWidget)`
            * :python:`wavelength is not None and not isinstance(wavelength, np.ndarray)`
            * :python:`wavelength_model is not None and not isinstance(wavelength_model, np.ndarray)`
            * :python:`wavelength_unit is not None and not isinstance(wavelength_unit, str)`
            * 
        '''
        
        if not isinstance(parent, QWidget):
            raise TypeError(f'Parent widget has type {type(parent)} but it should be of type QWidget.')
            
        if wavelength is not None and not isinstance(wavelength, np.ndarray):
            raise TypeError(f'wavelength range has type {type(wavelength)} but it should be None or numpy.ndarray.')
            
        if wavelength_model is not None and not isinstance(wavelength_model, np.ndarray):
            raise TypeError(f'model wavelength range has type {type(wavelength_model)} but it should be None or numpy.ndarray.')
        
        if wavelength_unit is not None and not isinstance(wavelength_unit, str):
            raise TypeError(f'wavelength unit has type {type(wavelength_unit)} but it should be None or str.')
        
        #: Parent widget
        self.__parent = parent

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
    
    #####################################################
    #        Updating properties of the spectrum        #
    #####################################################
    
    def update_spectrum(self, spectrum: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the spectrum.
        
        :param spectrum: 1D spectrum to show
        :type spectrum: numpy.ndarray
        '''
        
        # Update the y-data of the spectrum
        self.spectrum = spectrum
        
        # If first instantiation of the artist
        if self.__spec_artist is None:
            
            if self.wavelength is not None:
                
                self.__spec_artist = self.ax.plot(self.wavelength, self.spectrum, color='k', lw=1.5)
                
                # Update the layout of the figure
                self.ax.spines['bottom'].set_visible(True)
                self.ax.tick_params(bottom=True, labelbottom=True)
                self.ax.set_xlabel('Wavelength')
                
            else:
                self.__spec_artist = self.ax.plot(self.spectrum, color='k', lw=1.5)
                
            self.__spec_artist[0].set_drawstyle('steps')
                
            
        # Otherwise, just update the y-data of the artist
        else:
            
            self.__spec_artist[0].set_ydata(self.spectrum)
            self.ax.set_ylim([0.9*np.nanmin(self.spectrum), 1.1*np.nanmax(self.spectrum)])
        
        return
    
    def update_model_spectrum(self, spectrum: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Update the model spectrum.
        
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
    
class Dock_widget_spectrum(QDockWidget):
    
    def __init__(self, parent: QWidget, 
                 wavelength       : np.ndarray = None,
                 wavelength_model : np.ndarray = None
                ):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget
        :type parent: PyQt6.QtWidgets.QWidget
        
        Keyword arguments
        -----------------
        
        :param wavelength: wavelength for the data spectrum
        :type wavelength: numpy.ndarray or :python:`None`
        :param wavelength_model: wavelength for the model spectrum
        :type wavelength_model: numpy.ndarray or :python:`None`
        '''
        
        self.__parent      = parent
        
        super().__init__()
        
        # Canvas holding the spectrum
        self.__spec_canvas = Mpl_spectrum_canvas(self, 
                                                 wavelength       = wavelength, 
                                                 wavelength_model = wavelength_model
                                                )
        
        self.setWidget(self.spec_canvas)
        
        # Disable the possibility to close the dock
        self.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetMovable)
        
    @property
    def parent(self) -> QWidget:
        
        return self.__parent
        
    @property
    def spec_canvas(self) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Canvas holding the spectrum
        '''
        
        return self.__spec_canvas