#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
.. function:: pbxplore.chains_from_files(path_list)

   See :func:`pbxplore.structure.chains_from_files`

.. function:: pbxplore.chains_from_trajectory(trajectory, topology)

   See :func:`pbxplore.structure.chains_from_trajectory` 

.. function:: pbxplore.assign(dihedrals)

   See :func:`pbxplore.assign.assign`
"""

from .structure import chains_from_files, chains_from_trajectory
from .assign import assign
from . import PB
from . import io
from . import structure
from . import analysis
