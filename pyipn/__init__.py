# -*- coding: utf-8 -*-

"""Top-level package for pyipn."""


from .universe import Universe
from .fit import Fit
from .stan_models import get_stan_model, plot_stan_fit
from .io.package_utils import copy_template
from .possion_gen import pulse
import ligo.skymap.plot


__all__ = ["Universe", "copy_template"]

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
