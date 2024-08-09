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
        self.__image      = self.open_file(self.__args.image, self.__args.image_ext, ndim=2)  
        
        # Load data cube
        self.__cube       = self.open_file(self.__args.cube, self.__args.cube_ext, ndim=3)  
        
        # Load cube model
        self.__cube_model = self.open_file(self.__args.cube_model, self.__args.cube_model_ext, ndim=3)
        
        if all([i is None for i in (self.image, self.cube, self.cube_model)]):
            raise IOError('At least one file must be provided. Either a 2D map with the -i (--image) tag, or the 3D data cube with -c (--cube) tag, or the 3D cube model with -m (--model) tag.')
        
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
            self.bottom_dock = Dock_widget_spectrum(self)
            
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
    def cube_model(self) -> np.ndarray | None:
        r'''
        .. codeauthor:: Wilfried Mercier - LAM <wilfried.mercier@lam.fr>
        
        3D cube_model.
        '''
        
        return self.__cube_model
    
    #########################
    #       IO methods      #
    #########################
    
    @staticmethod
    def open_file(file: str | None, ext: int,
                  ndim : int | None = None
                 ) -> np.ndarray:
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
        
        Keyword arguments
        -----------------
        
        :param int ndim: number of dimensions the data should have. Possibilities are 2 or 3. If None, no check is performed
        
        :returns: data from the file
        :rtype: numpy.ndarray
        
        :raises: 
            * :python:`TypeError` if :python:`not isinstance(file, str)` or :python:`not isinstance(ext, int)` or :python:`not isinstance(ndim, int)`
            * :python:`IOError` if **file** is not a correct file name or symbolic link
            * :python:`ValueError` if the file extension is not recognized or `ndim not in [2, 3]` or `data.ndim != ndim`
        '''
        
        if file is None:
            return None
        
        if ndim is not None and not isinstance(ndim, int):
            raise TypeError(f'ndim has type {type(ndim)} but it should have type int.')
        
        if ndim is not None and ndim not in [2, 3]:
            raise ValueError(f'ndim is equal to {ndim} but it must be 2 or 3.')
        
        if not isinstance(file, str):
            raise TypeError(f'File has type {type(file)} but it must be of type str.')
            
        if not isinstance(ext, int):
            raise TypeError(f'Extension number has type {type(ext)} but it must be of type int.')
            
        if not opath.isfile(file):
            raise IOError(f'File {file} not found.')
            
        if (file_ext := opath.splitext(file)[1].lower()) in ['.fit', '.fits']:
            
            with fits.open(file) as hdul:    
                data = hdul[ext].data
        
        else:
            raise ValueError(f'File type with extension {file_ext} not supported. Supported types are: fits and fit.')
            
        if ndim is not None and data.ndim != ndim:
            raise ValueError(f'Data has dimensions {data.ndim} but it should have a dimension of {ndim}.')
            
        return data
        

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
    parser.add_argument('--cube_ext',         type=int, nargs='?', default=0, help='File extension for the 2D map (if applicable).')
    parser.add_argument('--cube_model_ext',   type=int, nargs='?', default=0, help='File extension for the 2D map (if applicable).')
    
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