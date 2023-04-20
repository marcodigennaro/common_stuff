# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:56:54 2021

@author: MDI0316
"""

import numpy as np
from typing import Any


def morse(r: float, e: float = 1., a: float = 1., r_eq: float = 1.) -> Any:
    """
    Morse potential
    Args:
        r: radius
        e: epsilon value
        a: alpha
        r_eq: equilibrium distance
    Returns:
        Value of Morse potential at given distance
    """
    return -e * (1 - (1 - np.exp(-a * (r - r_eq))) ** 2)


def lennard_jones(r: float, e: float = 1., s: float = 1.) -> Any:
    """
    Lennard Jones 6-12 potential
    Args:
        r: distance
        e: epsilon value
        s: sigma parameter: where the potential is null
    Returns:
        Value of Lennard Jones potential at distance r
    """
    return 4 * e * ((s / r) ** 12 - (s / r) ** 6)
