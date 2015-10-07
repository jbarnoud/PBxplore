#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Input/Output in files --- :mod:`pbxplore.io`
============================================

Deals with writing and reading of files in various formats.

Fasta
-----

.. autofunction:: read_fasta

.. autofunction:: read_several_fasta

.. autofunction:: write_fasta

Results of an assignation
-------------------------

.. autofunction:: write_phipsi

.. autofunction:: write_flat

Results af analyses
-------------------

.. autofunction:: write_count_matrix

.. autofunction:: write_neq
"""

from .fasta import read_fasta, read_several_fasta, write_fasta
from .write import write_phipsi, write_flat, write_count_matrix, write_neq
