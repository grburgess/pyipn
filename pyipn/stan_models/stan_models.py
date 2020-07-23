import cmdstanpy
import arviz
import ipyvolume as ipv
import numpy as np
import pkg_resources
import os

_available_models = ["rff_ipn.stan", "rff_omega_ipn.stan", "rff_bw_ipn.stan"]


def get_stan_model(stan_model):

    assert (
        stan_model in _available_models
    ), f"{stan_model} is not in {','.join(_available_models)}"

    stan_file =  pkg_resources.resource_filename(
                    "pyipn", os.path.join("stan_models", stan_model))

    
    model = cmdstanpy.CmdStanModel(
        stan_file=stan_file, cpp_options={"STAN_THREADS": "TRUE"}
    )

    return model


def plot_stan_fit(fit, universe, cmap="Set1", color="blue"):

    ar = arviz.from_cmdstanpy(fit)

    raw_xyz = np.array(ar.posterior.grb_xyz).reshape(-1, ar.posterior.grb_xyz.shape[-1])

    rad = universe.grb_radius

    # scatter = rad + np.random.normal(0, rad * 0.05, size=len(raw_xyz))

    universe.plot_all_annuli(cmap=cmap, lw=3, threeD=True)

    xyz = rad * raw_xyz

    ipv.scatter(
        xyz[:, 0], xyz[:, 1], xyz[:, 2], marker="sphere", color=color, size=0.7
    )
