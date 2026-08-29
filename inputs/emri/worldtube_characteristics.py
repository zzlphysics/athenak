#!/usr/bin/env python3
"""Linearized incoming-mode boundary projection for an EMRI MHD worldtube.

The eigensystem is formed numerically from the ideal special-relativistic MHD
conservative variables and fluxes in a local orthonormal frame.  This follows the
local-characteristic route for GRMHD: geometry supplies an orthonormal tetrad at the
boundary point, while the seven-wave RMHD projection is performed in that frame.

This module is a reference/preprocessing implementation, not yet the AthenaK device
boundary hook.  It establishes the state convention, incoming-wave sign, degeneracy
checks, and manufactured tests that the later C++ implementation must reproduce.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

import worldtube_flux_emf as worldtube


# Worldtube ideal-gas state after the ordinary AthenaK primitive variables:
# rho,u1,u2,u3,pgas,bcc1,bcc2,bcc3.  Passive scalars, when present, precede bcc1.
IDEAL_MHD_STATE_VARIABLES = (
    "rho",
    "u1",
    "u2",
    "u3",
    "pgas",
    "bcc1",
    "bcc2",
    "bcc3",
)
FACE_PRIMITIVE_VARIABLES = (
    "rho",
    "pgas",
    "u_normal",
    "u_tangent_u",
    "u_tangent_v",
    "b_tangent_u",
    "b_tangent_v",
)


def _finite_vector(values: object, size: int, name: str) -> np.ndarray:
    result = np.asarray(values)
    if result.shape != (size,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must be a finite vector with length {size}")
    return result


def _validate_adiabatic_index(adiabatic_index: float) -> float:
    result = float(adiabatic_index)
    if not math.isfinite(result) or result <= 1.0:
        raise ValueError("adiabatic_index must be finite and greater than one")
    return result


def srmhd_conserved_flux(
    face_primitive: object,
    normal_magnetic: float,
    adiabatic_index: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return seven conserved variables and their outward-normal SRMHD flux.

    ``face_primitive`` is ``rho,p,u_n,u_u,u_v,B_u,B_v``.  AthenaK's velocity
    primitive is the spatial four-velocity, so ``W=sqrt(1+u_i u^i)``.  The normal
    magnetic field is a CT constraint/parameter and is not a seventh wave variable.
    """

    primitive = _finite_vector(face_primitive, 7, "face_primitive")
    gamma = _validate_adiabatic_index(adiabatic_index)
    normal_field = np.asarray(normal_magnetic, dtype=primitive.dtype)
    if not np.isfinite(normal_field).all():
        raise ValueError("normal_magnetic must be finite")
    density, pressure = primitive[:2]
    if np.isrealobj(primitive) and (density <= 0.0 or pressure <= 0.0):
        raise ValueError("density and pressure must be positive")
    spatial_u = primitive[2:5]
    magnetic = np.asarray(
        (normal_field, primitive[5], primitive[6]), dtype=primitive.dtype
    )
    lorentz = np.sqrt(1.0 + np.dot(spatial_u, spatial_u))
    velocity = spatial_u / lorentz
    enthalpy_density = density + gamma * pressure / (gamma - 1.0)
    magnetic_u = np.dot(magnetic, spatial_u)
    magnetic_four_time = magnetic_u
    magnetic_four_space = (
        magnetic + magnetic_four_time * spatial_u
    ) / lorentz
    magnetic_squared = (
        np.dot(magnetic, magnetic) + magnetic_u * magnetic_u
    ) / (lorentz * lorentz)
    total_enthalpy = enthalpy_density + magnetic_squared
    total_pressure = pressure + 0.5 * magnetic_squared

    mass = density * lorentz
    momentum = (
        total_enthalpy * lorentz * spatial_u
        - magnetic_four_time * magnetic_four_space
    )
    energy = (
        total_enthalpy * lorentz * lorentz
        - total_pressure
        - magnetic_four_time * magnetic_four_time
    )
    tau = energy - mass
    conserved = np.asarray(
        (mass, momentum[0], momentum[1], momentum[2], tau, magnetic[1], magnetic[2])
    )

    momentum_flux = (
        total_enthalpy * spatial_u[0] * spatial_u
        - magnetic_four_space[0] * magnetic_four_space
    )
    momentum_flux[0] += total_pressure
    flux = np.asarray(
        (
            mass * velocity[0],
            momentum_flux[0],
            momentum_flux[1],
            momentum_flux[2],
            momentum[0] - mass * velocity[0],
            magnetic[1] * velocity[0] - magnetic[0] * velocity[1],
            magnetic[2] * velocity[0] - magnetic[0] * velocity[2],
        )
    )
    return conserved, flux


def _complex_step_jacobian(
    function: object, primitive: np.ndarray, step: float = 1.0e-30
) -> np.ndarray:
    result = np.empty((primitive.size, primitive.size), dtype=np.float64)
    for column in range(primitive.size):
        perturbed = primitive.astype(np.complex128)
        perturbed[column] += 1j * step
        values = np.asarray(function(perturbed))
        result[:, column] = np.imag(values) / step
    if not np.isfinite(result).all():
        raise RuntimeError("complex-step RMHD Jacobian contains a non-finite value")
    return result


@dataclass(frozen=True)
class CharacteristicBasis:
    """Primitive-space RMHD eigenbasis ordered by increasing outward speed."""

    speeds: np.ndarray
    right_eigenvectors: np.ndarray
    left_eigenvectors: np.ndarray
    condition_number: float
    jacobian_residual: float


@dataclass(frozen=True)
class LinearHLLDiagnostics:
    """Linear per-mode response of an HLL boundary flux at one reference state.

    ``flux_gain`` compares the HLL flux perturbation with the physical flux
    perturbation for the same eigenmode.  An outgoing perturbation is placed on the
    interior side and an incoming perturbation on the exterior side.  A gain of one is
    exact.  Stationary modes have undefined gain and are reported as ``nan``.

    This is a local reflection-risk proxy, not a measured wave-packet reflection
    coefficient: the latter also depends on reconstruction, resolution, and timestep.
    """

    speeds: np.ndarray
    flux_gain: np.ndarray
    gain_error: np.ndarray
    minimum_signal_speed: float
    maximum_signal_speed: float


def characteristic_basis(
    face_primitive: object,
    normal_magnetic: float,
    adiabatic_index: float,
    imaginary_tolerance: float = 2.0e-8,
    condition_limit: float = 1.0e10,
) -> CharacteristicBasis:
    """Numerically construct the seven-wave primitive RMHD eigensystem."""

    primitive = np.asarray(
        _finite_vector(face_primitive, 7, "face_primitive"), dtype=np.float64
    )
    if primitive[0] <= 0.0 or primitive[1] <= 0.0:
        raise ValueError("density and pressure must be positive")
    gamma = _validate_adiabatic_index(adiabatic_index)
    if not math.isfinite(condition_limit) or condition_limit <= 1.0:
        raise ValueError("condition_limit must be finite and greater than one")

    def conserved(values: np.ndarray) -> np.ndarray:
        return srmhd_conserved_flux(values, normal_magnetic, gamma)[0]

    def flux(values: np.ndarray) -> np.ndarray:
        return srmhd_conserved_flux(values, normal_magnetic, gamma)[1]

    conserved_jacobian = _complex_step_jacobian(conserved, primitive)
    flux_jacobian = _complex_step_jacobian(flux, primitive)
    try:
        primitive_jacobian = np.linalg.solve(conserved_jacobian, flux_jacobian)
    except np.linalg.LinAlgError as error:
        raise RuntimeError(
            "RMHD primitive-to-conserved Jacobian is singular"
        ) from error
    eigenvalues, eigenvectors = np.linalg.eig(primitive_jacobian)
    imaginary_scale = max(float(np.max(np.abs(eigenvalues.real))), 1.0)
    imaginary_max = max(
        float(np.max(np.abs(eigenvalues.imag))),
        float(np.max(np.abs(eigenvectors.imag))),
    )
    if imaginary_max > imaginary_tolerance * imaginary_scale:
        raise RuntimeError(
            "RMHD characteristic basis is not numerically real; the state may be "
            "degenerate or non-hyperbolic"
        )
    order = np.argsort(eigenvalues.real)
    speeds = eigenvalues.real[order]
    right = eigenvectors.real[:, order]
    # Fix arbitrary column norms/signs so cross-language tests remain stable.
    for column in range(right.shape[1]):
        norm = float(np.linalg.norm(right[:, column]))
        if not math.isfinite(norm) or norm == 0.0:
            raise RuntimeError("RMHD characteristic eigenvector has zero norm")
        right[:, column] /= norm
        pivot = int(np.argmax(np.abs(right[:, column])))
        if right[pivot, column] < 0.0:
            right[:, column] *= -1.0
    condition = float(np.linalg.cond(right))
    if not math.isfinite(condition) or condition > condition_limit:
        raise RuntimeError(
            f"RMHD characteristic basis is ill-conditioned ({condition:.6g}); "
            "use a degeneracy-safe analytic basis before replaying this state"
        )
    left = np.linalg.inv(right)
    residual = float(
        np.max(np.abs(primitive_jacobian @ right - right * speeds[None, :]))
    )
    if np.any(speeds < -1.0 - 2.0e-10) or np.any(speeds > 1.0 + 2.0e-10):
        raise RuntimeError("RMHD characteristic speed lies outside the light cone")
    return CharacteristicBasis(speeds, right, left, condition, residual)


def linear_hll_mode_gains(
    face_primitive: object,
    normal_magnetic: float,
    adiabatic_index: float,
    speed_tolerance: float = 1.0e-10,
) -> LinearHLLDiagnostics:
    """Return the linear HLL boundary-flux gain for all seven RMHD modes.

    The outward normal points from the inner simulation into the supplied worldtube
    state.  Positive-speed perturbations therefore live on the interior side of the
    Riemann problem, while negative-speed perturbations live on the exterior side.
    The calculation uses the extremal speeds of the same reference state, matching the
    local linear limit of a two-wave HLL solver.
    """

    if not math.isfinite(speed_tolerance) or speed_tolerance < 0.0:
        raise ValueError("speed_tolerance must be finite and nonnegative")
    basis = characteristic_basis(face_primitive, normal_magnetic, adiabatic_index)
    return linear_hll_gains_from_speeds(basis.speeds, speed_tolerance)


def linear_hll_gains_from_speeds(
    characteristic_speeds: object,
    speed_tolerance: float = 1.0e-10,
) -> LinearHLLDiagnostics:
    """Evaluate the linear two-wave HLL response from seven ordered speeds."""

    if not math.isfinite(speed_tolerance) or speed_tolerance < 0.0:
        raise ValueError("speed_tolerance must be finite and nonnegative")
    speeds = np.asarray(
        _finite_vector(characteristic_speeds, 7, "characteristic_speeds"),
        dtype=np.float64,
    )
    if np.any(np.diff(speeds) < 0.0):
        raise ValueError("characteristic_speeds must be ordered")
    minimum = min(0.0, float(speeds[0]))
    maximum = max(0.0, float(speeds[-1]))
    gains = np.full(7, np.nan, dtype=np.float64)
    moving = np.abs(speeds) > speed_tolerance
    if minimum >= 0.0:
        gains[moving & (speeds > 0.0)] = 1.0
    elif maximum <= 0.0:
        gains[moving & (speeds < 0.0)] = 1.0
    else:
        denominator = maximum - minimum
        outgoing = moving & (speeds > 0.0)
        incoming = moving & (speeds < 0.0)
        gains[outgoing] = (
            maximum
            * (speeds[outgoing] - minimum)
            / (denominator * speeds[outgoing])
        )
        gains[incoming] = (
            minimum
            * (maximum - speeds[incoming])
            / (denominator * speeds[incoming])
        )
    errors = np.abs(gains - 1.0)
    return LinearHLLDiagnostics(
        speeds=speeds,
        flux_gain=gains,
        gain_error=errors,
        minimum_signal_speed=minimum,
        maximum_signal_speed=maximum,
    )


def face_basis(face_name: str) -> np.ndarray:
    """Return rows containing the outward normal, local +u, and local +v axes."""

    if face_name not in worldtube.ORIENTATIONS:
        raise ValueError(f"unknown cubical face {face_name}")
    orientation = worldtube.ORIENTATIONS[face_name]
    result = np.zeros((3, 3), dtype=np.float64)
    result[0, orientation.normal_axis] = orientation.normal_sign
    result[1, orientation.u_axis] = orientation.u_sign
    result[2, orientation.v_axis] = orientation.v_sign
    if not np.allclose(np.linalg.det(result), 1.0, rtol=0.0, atol=0.0):
        raise RuntimeError("worldtube face basis is not right handed")
    return result


def state_to_face_primitive(
    state: object, face_name: str, normal_magnetic: float
) -> np.ndarray:
    """Rotate an eight-variable source-frame state into seven face primitives."""

    source = np.asarray(_finite_vector(state, 8, "MHD state"), dtype=np.float64)
    if source[0] <= 0.0 or source[4] <= 0.0:
        raise ValueError("MHD state density and pressure must be positive")
    basis = face_basis(face_name)
    face_velocity = basis @ source[1:4]
    face_magnetic = basis @ source[5:8]
    face_magnetic[0] = float(normal_magnetic)
    return np.asarray(
        (
            source[0],
            source[4],
            face_velocity[0],
            face_velocity[1],
            face_velocity[2],
            face_magnetic[1],
            face_magnetic[2],
        )
    )


def face_primitive_to_state(
    primitive: object, face_name: str, normal_magnetic: float
) -> np.ndarray:
    """Rotate seven face primitives back to the source-frame eight-variable state."""

    face = np.asarray(
        _finite_vector(primitive, 7, "face_primitive"), dtype=np.float64
    )
    basis = face_basis(face_name)
    source_velocity = basis.T @ face[2:5]
    source_magnetic = basis.T @ np.asarray(
        (normal_magnetic, face[5], face[6]), dtype=np.float64
    )
    return np.asarray(
        (
            face[0],
            source_velocity[0],
            source_velocity[1],
            source_velocity[2],
            face[1],
            source_magnetic[0],
            source_magnetic[1],
            source_magnetic[2],
        )
    )


@dataclass(frozen=True)
class ProjectionDiagnostics:
    speeds: np.ndarray
    incoming: np.ndarray
    amplitudes: np.ndarray
    applied_fraction: float
    eigenvector_condition_number: float
    eigensystem_residual: float


def project_incoming_characteristics(
    interior_primitive: object,
    exterior_primitive: object,
    normal_magnetic: float,
    adiabatic_index: float,
    speed_tolerance: float = 1.0e-10,
    positivity_floor: float = 1.0e-14,
) -> tuple[np.ndarray, ProjectionDiagnostics]:
    """Replace only modes propagating inward relative to the outward face normal."""

    interior = np.asarray(
        _finite_vector(interior_primitive, 7, "interior_primitive"), dtype=np.float64
    )
    exterior = np.asarray(
        _finite_vector(exterior_primitive, 7, "exterior_primitive"), dtype=np.float64
    )
    if np.any(interior[:2] <= 0.0) or np.any(exterior[:2] <= 0.0):
        raise ValueError("interior and exterior density/pressure must be positive")
    if not math.isfinite(speed_tolerance) or speed_tolerance < 0.0:
        raise ValueError("speed_tolerance must be finite and nonnegative")
    if not math.isfinite(positivity_floor) or positivity_floor <= 0.0:
        raise ValueError("positivity_floor must be finite and positive")

    reference = 0.5 * (interior + exterior)
    # A logarithmic mean for positive thermodynamic variables is symmetric and remains
    # physical even when the two sides differ by many orders of magnitude.
    reference[:2] = np.sqrt(interior[:2] * exterior[:2])
    basis = characteristic_basis(reference, normal_magnetic, adiabatic_index)
    difference = exterior - interior
    amplitudes = basis.left_eigenvectors @ difference
    incoming = basis.speeds < -speed_tolerance
    applied = basis.right_eigenvectors @ (incoming * amplitudes)

    fraction = 1.0
    for variable in (0, 1):
        if interior[variable] + applied[variable] < positivity_floor:
            if applied[variable] >= 0.0:
                raise RuntimeError("characteristic positivity limiter is inconsistent")
            permitted = (interior[variable] - positivity_floor) / (-applied[variable])
            fraction = min(fraction, max(0.0, 0.99 * permitted))
    boundary = interior + fraction * applied
    if not np.isfinite(boundary).all() or np.any(boundary[:2] <= 0.0):
        raise RuntimeError(
            "incoming characteristic projection produced an invalid state"
        )
    diagnostics = ProjectionDiagnostics(
        speeds=basis.speeds,
        incoming=incoming,
        amplitudes=amplitudes,
        applied_fraction=fraction,
        eigenvector_condition_number=basis.condition_number,
        eigensystem_residual=basis.jacobian_residual,
    )
    return boundary, diagnostics


def characteristic_boundary_state(
    interior_state: object,
    exterior_state: object,
    face_name: str,
    normal_magnetic: float,
    adiabatic_index: float,
    speed_tolerance: float = 1.0e-10,
) -> tuple[np.ndarray, ProjectionDiagnostics]:
    """Project source-frame ideal-MHD state data at one cubical boundary cell."""

    interior_face = state_to_face_primitive(
        interior_state, face_name, normal_magnetic
    )
    exterior_face = state_to_face_primitive(
        exterior_state, face_name, normal_magnetic
    )
    boundary_face, diagnostics = project_incoming_characteristics(
        interior_face,
        exterior_face,
        normal_magnetic,
        adiabatic_index,
        speed_tolerance=speed_tolerance,
    )
    return (
        face_primitive_to_state(boundary_face, face_name, normal_magnetic),
        diagnostics,
    )
