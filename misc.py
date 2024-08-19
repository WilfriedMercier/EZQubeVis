#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
..codeauthor:: Mercier Wilfried - LAM <wilfried.mercier@lam.fr>

Miscellaneous classes.
"""

import enum
import numpy as np

class Application_states(enum.Enum):
    
    LOCK = enum.auto()

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