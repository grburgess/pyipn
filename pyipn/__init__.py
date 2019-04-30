# -*- coding: utf-8 -*-

"""Top-level package for pyipn."""

__author__ = """J. Michael Burgess"""
__email__ = "jburgess@mpe.mpg.de"
__version__ = "0.1.0"


from .universe import Universe
from .io.package_utils import copy_template

__all__ = ["Universe", "copy_template"]
