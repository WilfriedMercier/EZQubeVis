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
from   PyQt6.QtCore                        import Qt
from   matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Custom import
from   misc                               import Application_states, DummyMouseEvent

class Mpl_im_canvas(FigureCanvas):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom matplotlib canvas that can hold an image.
    
    Heavily inspired by the mplCanvas class in PyQubeVis.
    '''
    
    def __init__(self, 
                 parent : QWidget,
                 root   : QWidget,
                 cmap   : str,
                ) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param parent: parent widget holding this widget
        :type parent: PyQt6.QtWidgets.QWidget
        :param root: root widget at the base of the application
        :type root: PyQt6.QtWidgets.QWidget
        :param str cmap: initial colormap used
        '''
        
        self.__parent = parent
        
        self.__root   = root
        
        # Using the setter defined below to automatically perform the checks
        #: Zoom strength. That's the multiplicative factor used to zoom in or zoom out with the scroll wheel.
        self.zoom_strength = 0.5
        
        # Mouse coordinates continuously read
        self.__mouse_coordinates = ()

        # Set the figure and the axis
        # Note that figure cannot be private because of the parent class that requires it to be public
        self.figure   = mplf.Figure()
        
        super().__init__(figure=self.figure)
        
        # Axis is set private and has no setter
        self.__ax     = self.figure.add_subplot(111)
        self.ax.set_axis_off()
        
        # Setup layout of the figure and axes
        self.figure.subplots_adjust(left=0, bottom=0., right=1, top=1)
        
        # im_artist is set private and has neither getter nor setter
        # Initialization to None before creating a first artist in self.update_image when an image is provided
        self.__im_artist = None
        
        # By default, we hide the cursor. Only in lock mode, do we see it
        self.setCursor(Qt.CursorShape.BlankCursor)
        
        ######################################
        #           Show the image           #
        ######################################
        
        # Based on what is given by the user, we show the image or, if not provided, the cube or, if not provided, the cube model
        if self.root.image is not None:
            image = self.root.image
        elif self.root.cube is not None:
            image = self.root.cube[0, :, :]
        elif self.root.cube_model is not None:
            image = self.root.cube_model[0, :, :]
        else:
            raise ValueError('No file provided, cannot show an image.')
            
        # Using the setter defined below to automatically perform the checks
        self.cmap     = cmap
        
        # Update the image
        self.update_image(image)
        self.draw()
        
        # Update the vmin and vmax values of the image
        if self.array is None:
            self.vmin = None
            self.vmax = None
        else:
            self.vmin = np.nanquantile(self.array, 0.1)
            self.vmax = np.nanquantile(self.array, 0.9)
        
        # Artist that adds a square around the currently selected pixel
        self.__highlight_rect = None
        
        ##################################
        #         Mask rectangle         #
        ##################################
        
        # The mask is equal to True for masked pixels and False otherwise
        
        # Artist that adds an expandable rectangle used to mask pixel values and a second artist used to show a filled area
        self.__mask_rect      = None
        self.__mask_rect_fill = None
        
        # Coordinates of the four corners of the mask used to update the mask on screen
        self.__mask_coord     = None
        
        # True and False mask: True for masked pixels
        self.__mask           = None if self.array is None else np.full(self.array.shape, False, dtype=bool)
        
        # Stack containing the coordinates of the pixels that were removed. Used to 
        self.__mask_stack     = []
        
        #######################################
        #           Events handling           #
        #######################################
        
        self.mpl_connect('button_press_event',   self.mouse_press)
        self.mpl_connect('button_release_event', self.mouse_release)
        self.mpl_connect('motion_notify_event',  self.mouse_move)
        self.mpl_connect('scroll_event',         self.scroll)
            
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
    def highlight_rect(self) -> matplotlib.lines.Line2D | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Line2D object used to highlight the currently selected pixel.
        '''
    
        return self.__highlight_rect
    
    @property
    def mask_rect(self) -> matplotlib.lines.Line2D | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Line2D object used to highlight the pixels that the user wants to mask. This corresponds to the edges of the mask zone.
        '''
    
        return self.__mask_rect
    
    @property
    def mask_rect_fill(self) -> matplotlib.collections.PolyCollection | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Artist object used to highlight the pixels that the user wants to mask. This corresponds to the filled area.
        '''
    
        return self.__mask_rect_fill
    
    @property
    def mask(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Pixel mask.
        '''
    
        return self.__mask
    
    @property
    def mask_coord(self) -> tuple[float] | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Coordinates of the four corner of the mask.
        '''
    
        return self.__mask_coord
    
    @property
    def mask_stack(self) -> list | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Pixel mask.
        '''
    
        return self.__mask_stack
    
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
        
        Set the cmap to a new one and update the image.
        
        .. note::
            
            This will not update the cmap on the image because the draw() method is not called here.
        
        :param str cmap: new colormap for the image
       
        :raises:
            * :python:`TypeError` if `not isinstance(cmap, str)`
            * :python:`ValueError` if the given cmap does not belong to the list of cmaps from matplotlib.pyplot
        '''
        
        if not isinstance(cmap, str):
            raise TypeError(f'cmap is of type {type(cmap)} but it must be of type str.')
        
        if cmap not in self.root.cmaps_ok:
            raise ValueError(f'cmap is {cmap} which does not belong to the following list of cmaps from matplotlib: {self.root.cmaps_ok}.')
        
        self.__cmap   = cmap
        
        # Update the cmap of the image
        if self.__im_artist is not None:
            self.__im_artist.set_cmap(self.cmap)

        return
    
    @property
    def mouse_coordinates(self) -> tuple[float]:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Position of the mouse cursor in figure coordinates.
        '''
        
        return self.__mouse_coordinates
    
    @property
    def zoom_strength(self) -> str:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Strength of the zoom used for zooming in and out.
        '''
        
        return self.__zoom_strength
    
    @zoom_strength.setter
    def zoom_strength(self, zoom_strength: int | float) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Set the zoom strength to a new value.
        
        :param zoom_strength: new strength of the zoom. Must be > 1
       
        :raises:
            * :python:`TypeError` if `not isinstance(cmap, str)`
            * :python:`ValueError` if the given cmap does not belong to the list of cmaps from matplotlib.pyplot
        '''
        
        if not isinstance(zoom_strength, (int, float)):
            raise TypeError(f'zoom_strength is of type {type(zoom_strength)} but it must be of type int or float.')
        
        if zoom_strength <= 0:
            raise ValueError(f'zoom_strength equal to {zoom_strength} but it must be strictly larger than 0.')
        
        self.__zoom_strength = zoom_strength

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
    
    def update_spectrum(self, xpos: float, ypos: float) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given x and y coordinates, update the spectrum.
        
        :param float xpos: x position to extract the spectrum from
        :param float ypos: y position to extract the spectrum from
        '''
        
        # Canvas holding the spectrum to be updated
        spec_canvas = self.root.bottom_dock.spec_canvas
        
        # Only do the spectrum update if a data cube is provided
        if self.root.cube is not None:
            
            # Update the data spectrum
            spec_canvas.update_spectrum(self.root.cube[:, ypos, xpos])
            
        # Only do the spectrum update if a model cube is provided
        if self.root.cube_model is not None:
            spec_canvas.update_model_spectrum(self.root.cube_model[:, ypos, xpos])
        
        return
    
    def move_highlight_rectangle(self, xpos, ypos) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given x and y coordinates, update the position of the highlight rectangle.
        
        :param float xpos: x position to place the highlight rectangle
        :param float ypos: y position to place the highlight rectangle
        '''
        
        xvals = [xpos - 0.5, xpos + 0.5, xpos + 0.5, xpos - 0.5, xpos - 0.5]
        yvals = [ypos - 0.5, ypos - 0.5, ypos + 0.5, ypos + 0.5, ypos - 0.5]
        
        if self.highlight_rect is None:
            self.__highlight_rect, = self.ax.plot(xvals, yvals, lw=3, color='k')
        else:
            self.highlight_rect.set_xdata(xvals)
            self.highlight_rect.set_ydata(yvals)
            
        return
    
    def remove_highlight_rectangle(self) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Remove the highlight rectangle from the figure.
        '''
        
        if self.__highlight_rect is not None:
            
            self.__highlight_rect.remove()
            self.__highlight_rect = None
        
        return
    
    def remove_mask_rectangle(self) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Remove the mask rectangle from the figure.
        '''
        
        if self.__mask_rect is not None:
            
            self.__mask_rect.remove()
            self.__mask_rect      = None
            
            
        if self.__mask_rect_fill is not None:
            
            self.__mask_rect_fill.remove()
            self.__mask_rect_fill = None
            
        return
        
    
    def update_mask_rectangle(self, xpos, ypos, just_move: bool = False) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given x and y coordinates, update the mask rectangle bounds.
        
        :param float xpos: x bound of the rectangle
        :param float ypos: y bound of the rectangle
        
        Keyword arguments
        -----------------
        
        :param bool just_move: whether to just move the rectangle or expand it
        '''
        
        if self.mask_rect is None or just_move:
            
            xvals = [xpos - 0.5, xpos + 0.5, xpos + 0.5, xpos - 0.5, xpos - 0.5]
            yvals = [ypos - 0.5, ypos - 0.5, ypos + 0.5, ypos + 0.5, ypos - 0.5]
            
            self.__mask_coord = (xpos - 0.5, ypos - 0.5)
            
            self.__mask_rect, = self.ax.plot(xvals, yvals, lw=3, color='darkblue')
            
        else:
            
            # If the filled part of the rectangle is not shown yet, we show it with dummy coordinates before updating them
            if self.__mask_rect_fill is not None:
                self.__mask_rect_fill.remove()
                
            if xpos >= self.mask_coord[0] + 0.5:
                self.mask_rect.set_xdata([init := self.mask_coord[0], end := xpos + 0.5, end, init, init])
            elif xpos <= self.mask_coord[0] + 0.5:
                self.mask_rect.set_xdata([init := self.mask_coord[0] + 1, end := xpos - 0.5, end, init, init])
            
            xcoord = [init, end]
                
            if ypos >= self.mask_coord[1] + 0.5:
                self.mask_rect.set_ydata([init := self.mask_coord[1], init, end := ypos + 0.5, end, init])
            elif ypos <= self.mask_coord[1] + 0.5:
                self.mask_rect.set_ydata([init := self.mask_coord[1] + 1, init, end := ypos - 0.5, end, init])
            
            self.__mask_rect_fill = self.ax.fill_between(xcoord, init, end)
            
            self.mask_rect_fill.set(color='royalblue', alpha=0.5, lw=0)
        return
    
    #########################################
    #           Mouse interaction           #
    #########################################
    
    def scroll(self, event) -> None:
        
        # Multiplicative factor
        step = self.zoom_strength
        
        # Invert scaling if zooming out
        if event.button == 'down':
            step *= -1
            
        xmin = self.ax.get_xlim()[0] + step
        xmax = self.ax.get_xlim()[1] - step
        
        ymin = self.ax.get_ylim()[0] + step
        ymax = self.ax.get_ylim()[1] - step
        
        if (xmax - xmin) < 1 or (ymax - ymin) < 1:
            return
    
        self.ax.set_xlim([xmin, xmax])
        self.ax.set_ylim([ymin, ymax])
        
        event.canvas.draw()
        
        return
    
    def mouse_press(self, event) -> None:
        
        # Handle left click
        if event.button == 1 and Application_states.MASK in self.root.states:
            
            # Activate mask-on state: masking can now start
            self.root.states.add(Application_states.MASK_ON)
            
            # Send a dummy event to update the shape of the mask rectangle
            event = DummyMouseEvent(self, self.__mouse_coordinates[0], self.__mouse_coordinates[1])
            self.mouse_move(event)
            
            # Send a statusbar message
            self.root.status_bar.showMessage('Masking on: keep the left click pressed and move the mouse or use keys to select pixels to mask.')
        
        return
    
    def mouse_release(self, event) -> None:
        
        # Handle left click
        if event.button == 1 and Application_states.MASK in self.root.states:
        
            # Deactivate mask-on state: masking is finished
            self.root.states.remove(Application_states.MASK_ON)
            
            # Send a dummy event to update the shape of the mask rectangle
            event = DummyMouseEvent(self, self.__mouse_coordinates[0], self.__mouse_coordinates[1])
            self.mouse_move(event)
            
            # Send a statusbar message
            self.root.status_bar.showMessage('Masking off: mask mode still activated. Press "m" to deactivate mask mode.')
    
        return
            
    
    def mouse_move(self, event) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Function called when the mouse is moved. It updates the spectrum at the current cursor location.
        '''
        
        # If None, we are outside and we do not update. If Lock state, we do nothing.
        if event.xdata is None or event.ydata is None:
            return
        
        # Store mouse coordinates
        self.__mouse_coordinates = (event.xdata, event.ydata)
        
        xpos = int(np.round(event.xdata))
        ypos = int(np.round(event.ydata))
        
        # If position is out of bounds, we do nothing
        if xpos < 0 or xpos >= self.root.cube.shape[2] or ypos < 0 or ypos >= self.root.cube.shape[1]:
            return
            
        ##########################################
        #            Rectangle updates           #
        ##########################################
        
        if Application_states.MASK in self.root.states:
            
            self.remove_highlight_rectangle()
            
            # If mask mode activated but the mask is not on yet, we only move the square around
            if Application_states.MASK_ON not in self.root.states:
        
                self.remove_mask_rectangle()
                self.update_mask_rectangle(xpos, ypos, just_move=True)
                
            # If mask is on, we update the mask as the mouse moves
            else:
                self.update_mask_rectangle(xpos, ypos, just_move=False)
            
        # If lock mode, moving the mouse does nothing
        elif Application_states.LOCK in self.root.states:
            
            return
        
        # Move the position of the highlight rectangle if not in lock state
        elif Application_states.LOCK not in self.root.states:
            
            self.remove_mask_rectangle()
            self.move_highlight_rectangle(xpos, ypos)
        
        #######################################
        #           Spectrum update           #
        #######################################
        
        # Update spectrum
        self.update_spectrum(xpos, ypos)
        
        # Apply changes to the spectrum
        self.root.bottom_dock.spec_canvas.draw()
        
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
                 parent           : QWidget,
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