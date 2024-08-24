#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Custom classes built around the FigureCanvas class from matplotlib.
"""

import matplotlib
import numpy                              as     np
import matplotlib.figure                  as     mplf
from   PyQt6.QtWidgets                    import QWidget, QTabWidget
from   PyQt6.QtCore                       import Qt
from   matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Custom import
from   .misc                              import (Application_states, Dummy_mouse_event, 
                                                  Base_widget_skeleton)

class Mpl_im_canvas(Base_widget_skeleton, FigureCanvas):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom matplotlib canvas that can hold an image.
    
    Heavily inspired by the mplCanvas class in PyQubeVis.
    '''
    
    def __init__(self, 
                 parent   : QWidget,
                 root     : QWidget,
                 data     : np.ndarray,
                 *args, **kwargs
                ) -> None:
        r'''    
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param data: object used to get the data from. Must be either a 2D or 3D array.
        :type data: numpy.ndarray
        '''

        # Note that figure cannot be private because of the parent class that requires it to be public
        self.figure   = mplf.Figure()
        
        Base_widget_skeleton.__init__(self, parent, root)
        FigureCanvas.__init__(self, *args, figure=self.figure, **kwargs)
        
        # Using the setter defined below to automatically perform the checks
        #: Zoom strength.
        self.zoom_strength = 0.5
        
        # Mouse coordinates continuously read
        self.__mouse_coordinates = ()
        
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
        
        # Update the image
        self.update_image(data)
        self.draw()
        
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
    
    def update_image(self, image: np.ndarray) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param image: image to show
        :type image: numpy.ndarray
        '''
        
        # Update the array, vmin and vmax values
        self.array = image
        
        # If first instantiation of the artist
        if self.im_artist is None:
            self.__im_artist = self.ax.imshow(self.array, 
                                              cmap   = self.parent.cmap, 
                                              vmin   = self.vmin, 
                                              vmax   = self.vmax,
                                              origin = 'lower'
                                             )
            
        # Otherwise, just update the data of the artist
        else:
            self.im_artist.set_data(self.array)
        
        return
    
    def update_spectrum_given_coordinates(self, xpos: float, ypos: float) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given x and y coordinates, update the spectrum.
        
        :param xpos: x position to extract the spectrum from
        :type xpos: :python:`float`
        :param ypos: y position to extract the spectrum from
        :type ypos: :python:`float`
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
    
    ###################################
    #       Highlight rectangle       #
    ###################################
    
    @property
    def highlight_rect(self) -> matplotlib.lines.Line2D | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Line2D object used to highlight the currently selected pixel.
        '''
    
        return self.__highlight_rect
    
    def move_highlight_rectangle(self, xpos: float, ypos: float) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given x and y coordinates, update the position of the highlight rectangle.
        
        :param xpos: x position to place the highlight rectangle
        :type xpos: :python:`float`
        :param ypos: y position to place the highlight rectangle
        :type ypos: :python:`float`
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
        
        if self.highlight_rect is not None:
            
            self.highlight_rect.remove()
            self.__highlight_rect = None
        
        return
    
    ##############################
    #       Mask rectangle       #
    ##############################
    
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
    
    def remove_mask_rectangle(self) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Remove the mask rectangle from the figure.
        '''
        
        # If the mask rectangle exists, we remove it from the figure and put it back to its default value of None
        if self.mask_rect is not None:
            
            self.mask_rect.remove()
            self.__mask_rect      = None
            
        # If the filled mask rectangle exists, we remove it from the figure and put it back to its default value of None
        if self.mask_rect_fill is not None:
            
            self.mask_rect_fill.remove()
            self.__mask_rect_fill = None
            
        return
        
    def update_mask_rectangle(self, xpos: float, ypos: float, just_move: bool = False) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given x and y coordinates, update the mask rectangle bounds.
        
        :param xpos: x bound of the rectangle
        :type xpos: :python:`float`
        :param ypos: y bound of the rectangle
        :type ypos: :python:`float`
        
        Keyword arguments
        -----------------
        
        :param just_move: whether to just move the rectangle or expand it
        :type just_move: :python:`bool`
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
    
    ###################
    #      Mask       #
    ###################
    
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
    
    ####################################################
    #       Data array and associated properties       #
    ####################################################
    
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
    
        :raises TypeError: if :python:`not isinstance(image, np.ndarray)`
        :raises ValueError: if :python:`image.ndim != 2`
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
        mini, maxi = float(np.nanquantile(self.array, 0.2)), float(np.nanquantile(self.array, 0.8))
        if mini >= maxi: 
            mini, maxi = float(np.nanmin(self.array)), float(np.nanmax(self.array))
        
        # Same value everywhere means we replace the array by NaN and min and max do not matter
        if mini >= maxi:    
            
            self.vmin    = -1.0
            self.vmax    =  1.0
            self.__array = image * np.nan
        
        self.vmin = mini
        self.vmax = maxi
        
        
        return
    
    @property
    def vmin(self) -> float:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Minimum value used to show the image.
        '''
        
        return self.__vmin
    
    @vmin.setter
    def vmin(self, value: float) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Set the minimum value used to show the image.
        
        :param value: value used for vmin
        :type value: :python:`float`
        
        :raises TypeError: if :python:`not isinstance(value, float)`
        '''
        
        if not isinstance(value, float):
            raise TypeError(f'Vmin value has type {type(value)} but it should be of type float.')
        
        self.__vmin = value
        
        return
    
    @property
    def vmax(self) -> float:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Maximum value used to show the image.
        '''
        
        return self.__vmax
    
    @vmax.setter
    def vmax(self, value: float) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Set the maximum value used to show the image.
        
        :param value: value used for vmax
        :type value: :python:`float`
        
        :raises TypeError: if :python:`not isinstance(value, float)`
        '''
        
        if not isinstance(value, float):
            raise TypeError(f'Vmax value has type {type(value)} but it should be of type float.')
        
        self.__vmax = value
        
        return
    
    #####################################
    #       Miscellaneous methods       #
    #####################################
    
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
    
    @property
    def im_artist(self) -> matplotlib.image.AxesImage | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Artist that contains the image shown.
        '''
        
        return self.__im_artist
    
    #########################################
    #           Mouse interaction           #
    #########################################
    
    def scroll(self, event) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Actions taken when the mouse is scrolled.
        '''
        
        # Multiplicative factor
        step = self.zoom_strength
        
        # Invert scaling if zooming out
        if event.button == 'down':
            step *= -1
            
        xmin = self.ax.get_xlim()[0] + step
        xmax = self.ax.get_xlim()[1] - step
        
        ymin = self.ax.get_ylim()[0] + step
        ymax = self.ax.get_ylim()[1] - step
        
        # Scrolling can invert when we reach the size of a pixel. In that case, we stop scrolling.
        if (xmax - xmin) < 1 or (ymax - ymin) < 1:
            return
    
        self.ax.set_xlim([xmin, xmax])
        self.ax.set_ylim([ymin, ymax])
        
        event.canvas.draw()
        
        return
    
    def mouse_press(self, event) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Actions taken when the mouse is pressed. This is used for masking.
        '''
        
        # Handle left click
        if event.button == 1 and Application_states.MASK in self.root.states:
            
            # Activate mask-on state: masking can now start
            self.root.states.add(Application_states.MASK_ON)
            
            # Send a dummy event to update the shape of the mask rectangle
            event = Dummy_mouse_event(self, self.__mouse_coordinates[0], self.__mouse_coordinates[1])
            self.mouse_move(event)
            
            # Send a statusbar message
            self.root.status_bar.showMessage('Masking on: keep the left click pressed and move the mouse or use keys to select pixels to mask.')
        
        return
    
    def mouse_release(self, event) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Actions taken when the mouse is realeased. This is used for masking.
        '''
        
        # Handle left click
        if event.button == 1 and Application_states.MASK in self.root.states:
            
            ####################################################################################
            #        Store positions of masked pixels and order of masking in the stack        #
            ####################################################################################
            
            xmin   = int(np.nanmin(self.mask_rect.get_data()[0]) + 0.5)
            xmax   = int(np.nanmax(self.mask_rect.get_data()[0]) + 0.5)
            xrange = range(xmin, xmax)
            
            ymin   = int(np.nanmin(self.mask_rect.get_data()[1]) + 0.5)
            ymax   = int(np.nanmax(self.mask_rect.get_data()[1]) + 0.5)
            yrange = range(ymin, ymax)
            
            tmp_stack = {}
            
            # Store (x, y) combinations and associated value in the array if the (x, y) combination does not already exist somewhere in the stack
            for x in xrange:
                for y in yrange:
                    for stack in self.mask_stack:
                        if (x, y) in stack.keys():
                            break
                    else:
                        tmp_stack[(x, y)] =  self.array[y, x]
                        
            if tmp_stack != {}:
                self.mask_stack.append(tmp_stack)
            
            self.array[[y for _ in xrange for y in yrange], [x for x in xrange for _ in yrange]] = np.nan
            self.update_image(self.array)
            self.draw()
        
            # Deactivate mask-on state: masking is finished
            self.root.states.remove(Application_states.MASK_ON)
            
            # Send a dummy event to update the shape of the mask rectangle
            event = Dummy_mouse_event(self, self.__mouse_coordinates[0], self.__mouse_coordinates[1])
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
        self.update_spectrum_given_coordinates(xpos, ypos)
        
        # Show the vertical line in the spectrum canvas if applicable (cube shown in the im canvas)
        self.root.bottom_dock.spec_canvas.update_spec_line_position(self.root.cube_pos)
        
        # Apply changes to the spectrum
        self.root.bottom_dock.spec_canvas.draw()
        
        # Apply changes to the image
        event.canvas.draw()
        
        return
    
class Tab_mpl_images(Base_widget_skeleton, QTabWidget):
    r'''
    .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
    
    Custom Tab widget that holds multiple matplotlib canvas.
    '''
    
    def __init__(self, parent: QWidget, root: QWidget, cmap: str, *args, **kwargs) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param cmap: colormap for the figure
        :type cmap: :python:`str`
        '''
        
        Base_widget_skeleton.__init__(self, parent, root)
        QTabWidget.__init__(self, *args, **kwargs)
        
        # Initialization to None
        self.__image_widget      = None
        self.__cube_widget       = None
        self.__cube_model_widget = None
        
        self.cmap                = cmap
        
        # Create one tab for the input image (if provided)
        if self.root.image is not None:
            
            self.__image_widget = Mpl_im_canvas(self, 
                                                self.root, 
                                                self.root.image
                                               )
            
            self.addTab(self.__image_widget, '2D &image')
            
        # Create one tab for the data cube (if provided)
        if self.root.cube is not None:
            
            self.__cube_widget       = Mpl_im_canvas(self, 
                                                     self.root, 
                                                     self.root.cube[self.root.cube_pos, :, :]
                                                    )
            
            self.addTab(self.__cube_widget, '&Data cube')
            
        # Create one tab for the cube model (if provided)
        if self.root.cube_model is not None:
            
            self.__cube_model_widget = Mpl_im_canvas(self, 
                                                     self.root, 
                                                     self.root.cube_model[self.root.cube_pos, :, :]
                                                    )
            
            self.addTab(self.__cube_model_widget, 'Cube &model')
            
        # Whatever what page was intialized we set it to the first one
        self.setCurrentIndex(0)
        
        return
    
    @property
    def image_widget(self) -> Mpl_im_canvas:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Widget holding the 2D image.
        '''
        
        return self.__image_widget
    
    @property
    def cube_widget(self) -> Mpl_im_canvas:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Widget holding the data cube.
        '''
        
        return self.__cube_widget
    
    @property
    def cube_model_widget(self) -> Mpl_im_canvas:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Widget holding the cube model.
        '''
        
        return self.__cube_model_widget
    
    ########################
    #       Colormap       #
    ########################
    
    @property
    def cmap(self) -> str:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Current colormap used to show the figures. The cmap is shared across all figures.
        '''
        
        return self.__cmap
    
    @cmap.setter
    def cmap(self, cmap: str) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Set the cmap to a new one and update the image.
        
        .. note::
            
            This will not update the cmap on the image because the draw() method is not called here.
        
        :param cmap: new colormap for the image
        :type cmap: :python:`str`
       
        :raises TypeError: if :python:`not isinstance(cmap, str)`
        :raises ValueError: if the given cmap does not belong to the allowed list of cmaps from matplotlib.pyplot
        '''
        
        if not isinstance(cmap, str):
            raise TypeError(f'cmap is of type {type(cmap)} but it must be of type str.')
        
        if cmap not in self.root.cmaps_ok:
            raise ValueError(f'cmap is {cmap} which does not belong to the following list of cmaps from matplotlib: {self.root.cmaps_ok}.')
        
        self.__cmap   = cmap
        
        # Update the cmap of the image
        if self.image_widget is not None:
            self.image_widget.im_artist.set_cmap(self.cmap)
            
        # Update the cmap of the cube
        if self.cube_widget is not None:
            self.cube_widget.im_artist.set_cmap(self.cmap)
            
        # Update the cmap of the cube model
        if self.cube_model_widget is not None:
            self.cube_model_widget.im_artist.set_cmap(self.cmap)

        return