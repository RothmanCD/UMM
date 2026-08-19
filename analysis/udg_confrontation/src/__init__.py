"""UMM first quantitative confrontation with DM-deficient UDG kinematics."""

from .physics import (
    A_FROZEN,
    A_FROZEN_UNC,
    N_G,
    SIGMA_REF,
    ALPHA_EFF,
    a_ind_powerlaw,
    a_ind_logistic,
    mean_sigma_within_re,
    wolf_sigma_from_mass,
    wolf_mass_from_sigma,
    a_bar_half_light,
    predict_sigma_three_models,
    G_KPC,
)

__all__ = [
    "A_FROZEN",
    "A_FROZEN_UNC",
    "N_G",
    "SIGMA_REF",
    "ALPHA_EFF",
    "a_ind_powerlaw",
    "a_ind_logistic",
    "mean_sigma_within_re",
    "wolf_sigma_from_mass",
    "wolf_mass_from_sigma",
    "a_bar_half_light",
    "predict_sigma_three_models",
    "G_KPC",
]
