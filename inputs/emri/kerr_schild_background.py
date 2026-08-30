"""Aligned-spin Cartesian Kerr-Schild background helpers.

This module is intentionally independent of the worldtube and frame builders so
both numerical-ADM and state-only analytic-Kerr source adapters can use one
implementation of the primary metric.
"""

from __future__ import annotations

import math

import numpy as np


def covariant_metric(
    position: object, primary_mass: float, primary_spin: float
) -> np.ndarray:
    """Return the unregularized Cartesian Kerr-Schild covariant metric."""

    point = np.asarray(position, dtype=np.float64)
    if point.shape != (3,) or not np.isfinite(point).all():
        raise ValueError("Kerr-Schild position must contain three finite values")
    spin2 = primary_spin**2
    radius2 = float(point @ point)
    adotx = primary_spin * point[2]
    radial_term = radius2 - spin2
    kerr_r2 = 0.5 * (
        radial_term + math.sqrt(radial_term**2 + 4.0 * adotx**2)
    )
    if primary_mass <= 0.0 or kerr_r2 <= 0.0:
        raise ValueError("Kerr-Schild metric point is singular or mass is invalid")
    kerr_r = math.sqrt(kerr_r2)
    denominator = kerr_r2**2 + adotx**2
    hfun = primary_mass * kerr_r2 * kerr_r / denominator
    spatial_denominator = kerr_r2 + spin2
    x_cross_a = np.asarray(
        (primary_spin * point[1], -primary_spin * point[0], 0.0)
    )
    null = np.empty(4, dtype=np.float64)
    null[0] = 1.0
    null[1:] = (
        kerr_r * point
        + x_cross_a
        + (adotx / kerr_r) * np.asarray((0.0, 0.0, primary_spin))
    ) / spatial_denominator
    metric = np.diag((-1.0, 1.0, 1.0, 1.0))
    metric += 2.0 * hfun * np.outer(null, null)
    return metric


def adm_fields(
    positions: object, primary_mass: float, dimensionless_spin: float
) -> np.ndarray:
    """Return analytic ADM fields in the standard ten-field worldtube order."""

    points = np.asarray(positions, dtype=np.float64)
    if (
        points.ndim != 2
        or points.shape[1] != 3
        or not np.isfinite(points).all()
        or not math.isfinite(primary_mass)
        or primary_mass <= 0.0
        or not math.isfinite(dimensionless_spin)
        or abs(dimensionless_spin) > 1.0
    ):
        raise ValueError("analytic Kerr ADM sample parameters are invalid")
    result = np.empty((points.shape[0], 10), dtype=np.float64)
    primary_spin = primary_mass * dimensionless_spin
    for index, point in enumerate(points):
        metric = covariant_metric(point, primary_mass, primary_spin)
        gamma = metric[1:, 1:]
        beta_lower = metric[0, 1:]
        beta = np.linalg.solve(gamma, beta_lower)
        lapse_squared = float(beta_lower @ beta - metric[0, 0])
        if lapse_squared <= 0.0 or np.min(np.linalg.eigvalsh(gamma)) <= 0.0:
            raise RuntimeError("analytic Kerr ADM decomposition is invalid")
        result[index] = (
            gamma[0, 0],
            gamma[0, 1],
            gamma[0, 2],
            gamma[1, 1],
            gamma[1, 2],
            gamma[2, 2],
            math.sqrt(lapse_squared),
            beta[0],
            beta[1],
            beta[2],
        )
    return result
