# -*- coding: utf-8 -*-

"""Top-level package for pyipn."""


from .universe import Universe
from .stan_models import get_stan_model, plot_stan_fit
from .io.package_utils import copy_template


__all__ = ["Universe", "copy_template"]

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
