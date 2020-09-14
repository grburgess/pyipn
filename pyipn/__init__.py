# -*- coding: utf-8 -*-

"""Top-level package for pyipn."""


from .universe import Universe
from .fit import Fit
from .stan_models import get_stan_model, list_stan_models
from .lightcurve import BinnedLightCurve

from .correlation import Correlator
from .io.package_utils import copy_template
from .possion_gen import pulse
import ligo.skymap.plot


__all__ = ["Universe", "copy_template", "Fit", "get_stan_model", "list_stan_models"]

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
