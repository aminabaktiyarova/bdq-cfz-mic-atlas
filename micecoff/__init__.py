"""
Interval-censored MIC estimation and epidemiological cut-off derivation.

A broth microdilution reading is an interval of concentrations. This package
fits the distribution those intervals came from, reports the shift between two
groups in doublings with an interval that resamples clusters of isolates, and
derives the cut-off above which a fitted wild-type population has negligible
mass.
"""

from .core import (  # noqa: F401
    check_doubling_series,
    ecoff,
    fit_censored_linear,
    fit_censored_normal,
    log_density_ratios,
    log_interval_mass,
    mic_bounds,
    mic_recorded,
    named_generator,
    parse_concentration,
    place_mics,
    resample_clusters,
    resample_shift,
    round_up_to_series,
    shift,
    simulate_reports,
    withhold_reason,
)

__version__ = "0.1.0"
