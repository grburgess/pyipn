import cmdstanpy
import arviz


_available_models = ["rff_ipn.stan", "rff_omega_ipn.stan"]


def get_stan_model(stan_model):

    assert (
        stan_model in _available_models
    ), f"{stan_model} is not in {','.join(_available_models)}"

    model = cmdstanpy.CmdStanModel(
        stan_file=stan_model, cpp_options={"STAN_THREADS": "TRUE"}
    )


    return model




