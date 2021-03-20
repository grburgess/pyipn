import cmdstanpy
import arviz
import ipyvolume as ipv
import numpy as np
import pkg_resources
import os

_available_models = [
    "rff_ipn.stan",
    "rff_omega_ipn.stan",
    "rff_omega_ipn_cos.stan",
    "rff_omega_ipn_earth.stan",
    "rff_bw_ipn.stan",
    "rff.stan",
    "rff_bw.stan",
    "rff_omega_ipn_dt.stan",
    "rff_omega.stan",
]


def get_stan_model(stan_model, mpi=False, threads=True):

    assert (
        stan_model in _available_models
    ), f"{stan_model} is not in {','.join(_available_models)}"

    stan_file = pkg_resources.resource_filename(
        "pyipn", os.path.join("stan_models", stan_model)
    )

    cpp_options = {}

    if mpi:
        cpp_options["STAN_MPI"] = True

    if threads:

        cpp_options["STAN_THREADS"] = True

    model = cmdstanpy.CmdStanModel(stan_file=stan_file, cpp_options=cpp_options)

    return model


def list_stan_models():
    for m in _available_models:

        print(m)
