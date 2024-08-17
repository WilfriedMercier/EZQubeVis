#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Inspired by PyQubeVis (https://gitlab.lam.fr/bepinat/PyQubeVis/-/tree/master?ref_type=heads) from Epinat Benoit (LAM).
"""

import sys
import signal
import argparse
import os.path                            as     opath
import numpy                              as     np
from   astropy.io                         import fits
from   PyQt6.QtWidgets                    import QApplication, QMainWindow
from   PyQt6.QtCore                       import Qt

# Custom imports
from   mpl_custom                        import Mpl_im_canvas, Dock_widget_spectrum

class Window(QMainWindow):
    
    def __init__(self, args: argparse.Namespace) -> None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        :param args: argument namespace built with the argparse library. Contains all the arguments that the program may use
        :type args: argparse.Namespace
        
        :raises:
            * :python:`TypeError` if :python:`not isinstance(args, argparse.Namespace)`
            * :python:`IOError` if no files were provided
        '''
        
        super().__init__()
        
        if not isinstance(args, argparse.Namespace):
            raise TypeError(f'args in Window instance has type {type(args)} but only argparse.Namespace is allowed.')
        
        #: Command line arguments
        self.__args   = args
        
        ###########################
        #        Load data        #
        ###########################
        
        # Load 2D map
        self.__image, self.__image_hdr = self.open_file(self.__args.image, 
                                                        self.__args.image_ext, 
                                                        self.__args.image_hdr_ext, 
                                                        ndim=2
                                                       )
        
        # Load data cube
        self.__cube, self.__cube_hdr = self.open_file(self.__args.cube, 
                                                      self.__args.cube_ext, 
                                                      self.__args.cube_hdr_ext,
                                                      ndim=3
                                                     )
        
        # Get spectral range of data cube
        self.__cube_wv    = self.compute_wavelength_range(self.__cube_hdr)
        
        # Load cube model
        self.__cube_model, self.__cube_model_hdr = self.open_file(self.__args.cube_model, 
                                                                  self.__args.cube_model_ext, 
                                                                  self.__args.cube_model_hdr_ext,
                                                                  ndim=3
                                                                 )
        
        # Get spectral range of cube model
        self.__cube_model_wv = self.compute_wavelength_range(self.__cube_model_hdr)
        
        # Check that there is at least one file provided
        if all([i is None for i in (self.image, self.cube, self.cube_model)]):
            raise IOError('At least one file must be provided. Either a 2D map with the -i (--image) tag, or the 3D data cube with -c (--cube) tag, or the 3D cube model with -m (--model) tag.')
        
        # Check that, if the wavelength ranges are computed, they match the number of spaxels in the cubes
        if self.__cube_wv is not None and self.__cube_wv.shape[0] != self.cube.shape[0]:
            print(f'WARNING: wavelength range derived from header has {self.cube_wv.shape[0]} spaxels but data cube has {self.cube.shape[0]} spaxels. This is inconsistent. No wavelength range will be shown.')
            self.__cube_wv = None
        
        if self.__cube_model_wv is not None and self.__cube_model_wv.shape[0] != self.cube_model.shape[0]:
            print(f'WARNING: wavelength range derived from header has {self.__cube_model_wv.shape[0]} spaxels but cube model has {self.__cube_model.shape[0]} spaxels. This is inconsistent. No wavelength range will be shown.')
            self.__cube_model_wv = None
            
        if self.cube is not None and self.cube_model is not None and (self.__cube_wv is None or self.__cube_model_wv is None):
            print('WARNING: Both data cube and cube model were provided but one of the wavelength ranges could not be properly computed. No wavelength range will be shown.')
            
            self.__cube_wv       = None
            self.__cube_model_wv = None
            
        ##################################################################
        #            Main frame that serves as central widget            #
        ##################################################################
        
        self.mpl_im_widget = Mpl_im_canvas(self, 'rainbow')
        
        self.setCentralWidget(self.mpl_im_widget)
        
        ###############################
        #            Docks            #
        ###############################
        
        # We only show the dock with the spectrum if 2 or 3 files are provided
        if self.cube is not None or self.cube_model is not None:
            
            # Add a dock at the bottom that will hold the spectrum
            self.bottom_dock = Dock_widget_spectrum(self, self.__cube_wv, self.__cube_model_wv)
            
            self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, 
                               self.bottom_dock, 
                               Qt.Orientation.Horizontal
                              )
            
            self.resizeDocks((self.bottom_dock, ), (300,), Qt.Orientation.Vertical)
        
        return
    
    ###################################
    #       Getters and setters       #
    ###################################
    
    @property
    def image(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        2D map.
        '''
        
        return self.__image
    
    @property
    def cube(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        3D data cube.
        '''
        
        return self.__cube
    
    @property
    def cube_hdr(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        3D data cube header.
        '''
        
        return self.__cube_hdr
    
    @property
    def cube_model(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        3D cube_model.
        '''
        
        return self.__cube_model
    
    @property
    def cube_model_hdr(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        3D cube model header.
        '''
        
        return self.__cube_model_hdr
    
    #########################
    #       IO methods      #
    #########################
    
    @staticmethod
    def compute_wavelength_range(hdr: fits.Header | None) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Given a header, compute the associated wavelength range if possible. If not possible, the spectrum will be shown without an associated wavelength range.
        
        :param hdr: header to extract the wavelength information from
        :type hdr: astropy.io.fits.Header
        
        :returns: wavelength range or None
        :rtype: numpy.ndarray or :python:`None`
        '''
        
        if hdr is None:
            return None
        
        if (key := 'CRVAL3') not in hdr or (key := 'CRPIX3') not in hdr or (key := 'NAXIS3') not in hdr:
            print(f'WARNING: {key} not in header. Cannot compute wavelength range.')
            return None
        
        return hdr['CRVAL3'] + hdr['CD3_3'] * np.linspace(init := hdr['CRPIX3'], init + (length := hdr['NAXIS3']), length)
        
    
    @staticmethod
    def open_file(file: str | None, ext: int, ext_hdr : int,
                  ndim : int | None = None
                 ) -> (np.ndarray, fits.Header | None):
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        Open a file.
        
        .. note::
            
            Currently supported files are:
                
                * FITS
                * FIT
        
        :param file: file to open. If None, no action is performed.
        :type file: :python:`str` or :python:`None`
        :param int ext: extension from which to get the file (e.g. for a FITS file). Not used if not applicable.
        :param int ext_hdr: extension from which to get the header information (e.g. for a FITS file). Not used if not applicable.
        
        Keyword arguments
        -----------------
        
        :param int ndim: number of dimensions the data should have. Possibilities are 2 or 3. If None, no check is performed
        
        :returns: data from the file and header (or None if not applicable)
        :rtype: (numpy.ndarray, astropy.io.fits.Header | :python:`None`)
        
        :raises: 
            * :python:`TypeError` if :python:`not isinstance(file, str)` or :python:`not isinstance(ext, int)` or :python:`not isinstance(ext_hdr, int)` or :python:`not isinstance(ndim, int)`
            * :python:`IOError` if **file** is not a correct file name or symbolic link
            * :python:`ValueError` if the file extension is not recognized or `ndim not in [2, 3]` or `data.ndim != ndim`
        '''
        
        if file is None:
            return None, None
        
        if ndim is not None and not isinstance(ndim, int):
            raise TypeError(f'ndim has type {type(ndim)} but it should have type int.')
        
        if ndim is not None and ndim not in [2, 3]:
            raise ValueError(f'ndim is equal to {ndim} but it must be 2 or 3.')
        
        if not isinstance(file, str):
            raise TypeError(f'File has type {type(file)} but it must be of type str.')
            
        if not isinstance(ext, int):
            raise TypeError(f'Extension number has type {type(ext)} but it must be of type int.')
            
        if not isinstance(ext_hdr, int):
            raise TypeError(f'Header extension number has type {type(ext_hdr)} but it must be of type int.')
            
        if not opath.isfile(file):
            raise IOError(f'File {file} not found.')
            
        if (file_ext := opath.splitext(file)[1].lower()) in ['.fit', '.fits']:
            
            with fits.open(file) as hdul:    
                data = hdul[ext].data
                hdr  = hdul[ext_hdr].header
        
        else:
            raise ValueError(f'File type with extension {file_ext} not supported. Supported types are: fits and fit.')
            
        if ndim is not None and data.ndim != ndim:
            raise ValueError(f'Data has dimensions {data.ndim} but it should have a dimension of {ndim}.')
            
        return data, hdr
        

def main(argv):
    
    ##################################################
    #             Command line arguments             #
    ##################################################
    
    # Instantiate the argument parser
    parser = argparse.ArgumentParser(prog        = 'EZQubeVis',
                                     description = 'EZQubeVis is a light-weight cube and 2D map vizualization tool inspired by PyQubeVis.'
                                    )
    
    # Mandatory arguments
    parser.add_argument('-i', '--image',      type=str, nargs='?', help='2D map to open.')
    parser.add_argument('-c', '--cube',       type=str, nargs='?', help='data cube to open')
    parser.add_argument('-m', '--cube_model', type=str, nargs='?', help='cube model to open')
    
    parser.add_argument('-s', '--save',       type=str, nargs='?', help='File within which the cleaned map will be saved when the program terminates.')
    
    parser.add_argument('--image_ext',        type=int, nargs='?', default=0, help='File extension for the 2D map (if applicable).')
    parser.add_argument('--cube_ext',         type=int, nargs='?', default=0, help='File extension for the 3D cube (if applicable).')
    parser.add_argument('--cube_model_ext',   type=int, nargs='?', default=0, help='File extension for the 3D cube model (if applicable).')
    
    parser.add_argument('--image_hdr_ext',      type=int, nargs='?', default=0, help='File extension for the header of the 2D map (if applicable).')
    parser.add_argument('--cube_hdr_ext',       type=int, nargs='?', default=0, help='File extension for the header of the 3D cube (if applicable).')
    parser.add_argument('--cube_model_hdr_ext', type=int, nargs='?', default=0, help='File extension for the header of the 3D cube model (if applicable).')
    
    # Get command-line arguments
    args   = parser.parse_args()

    #############################################
    #           Setting up the window           #
    #############################################
    
    # Main application and window
    app    = QApplication(sys.argv)
    
    # Handling SIGINT signal from terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    main_window = Window(args)
    main_window.show()
    
    # Add the possibility to close the window through the terminal
    sys.exit(app.exec())

if __name__ == '__main__':
    main(sys.argv[1:])