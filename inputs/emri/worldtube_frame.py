#!/usr/bin/env python3
"""Transform and constrain moving EMRI worldtube samples.

The magnetic data are treated as differential forms.  A spacetime map

    x^mu(T, X^a) = z^mu(T) + e^mu_a(T) X^a

pulls the Faraday two-form back with its full spacetime Jacobian.  In a
Cartesian chart with ``t=T`` this gives the edge electric one-form
``e_a dot (E + w cross B)``.  The motional term is therefore part of the
coordinate transformation, not an optional correction.

Interpolation and quadrature on a cut moving cube will not, in general, make
sampled endpoint fluxes and sampled edge EMFs obey the same discrete Faraday
law.  ``project_moving_samples`` minimally corrects the unique edge cochain so
that it has the requested curl on every surface face cell.  It never changes
the endpoint magnetic fluxes; a changing net flux through the closed cube is
rejected because no edge correction can repair it.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np

import worldtube_flux_emf as worldtube


FRAME_CONTRACT_CLASSIFICATION = "athenak-emri-affine-worldtube-frame-v1"


def _finite_array(values: object, shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite array with shape {shape}")
    return array


@dataclass(frozen=True)
class AffineFrame:
    """Instantaneous affine source-frame map and its time derivative.

    The columns of ``spatial_legs`` are ``e^mu_a``.  ``worldline_tangent`` is
    ``dz^mu/dT`` and ``spatial_leg_derivative`` is ``de^mu_a/dT``.  The time
    column of the spacetime Jacobian at local position ``X`` is consequently
    ``dz/dT + de_a/dT X^a``.
    """

    worldline_tangent: np.ndarray
    spatial_legs: np.ndarray
    spatial_leg_derivative: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "worldline_tangent",
            _finite_array(self.worldline_tangent, (4,), "worldline_tangent"),
        )
        object.__setattr__(
            self,
            "spatial_legs",
            _finite_array(self.spatial_legs, (4, 3), "spatial_legs"),
        )
        object.__setattr__(
            self,
            "spatial_leg_derivative",
            _finite_array(
                self.spatial_leg_derivative,
                (4, 3),
                "spatial_leg_derivative",
            ),
        )

    def jacobian(self, local_position: object = (0.0, 0.0, 0.0)) -> np.ndarray:
        """Return ``partial x^mu / partial X^A`` at one local position."""

        position = _finite_array(local_position, (3,), "local_position")
        result = np.empty((4, 4), dtype=np.float64)
        result[:, 0] = self.worldline_tangent + self.spatial_leg_derivative @ position
        result[:, 1:] = self.spatial_legs
        determinant = float(np.linalg.det(result))
        scale = max(float(np.linalg.norm(result, ord=2)), 1.0)
        if abs(determinant) <= 64.0 * np.finfo(float).eps * scale**4:
            raise ValueError("affine frame has a singular spacetime Jacobian")
        return result


def faraday_tensor(electric: object, magnetic: object) -> np.ndarray:
    """Build ``F_0i=-E_i, F_ij=epsilon_ijk B^k`` in a Cartesian chart."""

    electric_array = _finite_array(electric, (3,), "electric")
    magnetic_array = _finite_array(magnetic, (3,), "magnetic")
    result = np.zeros((4, 4), dtype=np.float64)
    result[0, 1:] = -electric_array
    result[1:, 0] = electric_array
    result[1, 2] = magnetic_array[2]
    result[2, 1] = -magnetic_array[2]
    result[2, 3] = magnetic_array[0]
    result[3, 2] = -magnetic_array[0]
    result[3, 1] = magnetic_array[1]
    result[1, 3] = -magnetic_array[1]
    return result


def electric_magnetic_from_faraday(tensor: object) -> tuple[np.ndarray, np.ndarray]:
    """Recover coordinate electric one-form and magnetic flux-density vector."""

    array = _finite_array(tensor, (4, 4), "faraday tensor")
    if not np.allclose(array + array.T, 0.0, rtol=0.0, atol=2.0e-13):
        raise ValueError("faraday tensor is not antisymmetric")
    electric = -array[0, 1:].copy()
    magnetic = np.asarray((array[2, 3], array[3, 1], array[1, 2]))
    return electric, magnetic


def pullback_two_form(tensor: object, jacobian: object) -> np.ndarray:
    """Pull a covariant two-form back through ``x=x(X)``."""

    form = _finite_array(tensor, (4, 4), "two-form")
    mapping = _finite_array(jacobian, (4, 4), "spacetime jacobian")
    return mapping.T @ form @ mapping


def pullback_electric_magnetic(
    electric: object,
    magnetic: object,
    frame: AffineFrame,
    local_position: object = (0.0, 0.0, 0.0),
) -> tuple[np.ndarray, np.ndarray]:
    """Return EMF one-form and magnetic flux density in affine coordinates."""

    pulled = pullback_two_form(
        faraday_tensor(electric, magnetic), frame.jacobian(local_position)
    )
    return electric_magnetic_from_faraday(pulled)


def transform_contravariant_vector(vector: object, jacobian: object) -> np.ndarray:
    """Transform ``V^mu`` to affine-coordinate components ``V^A``."""

    source = _finite_array(vector, (4,), "contravariant vector")
    mapping = _finite_array(jacobian, (4, 4), "spacetime jacobian")
    try:
        result = np.linalg.solve(mapping, source)
    except np.linalg.LinAlgError as error:
        raise ValueError("spacetime jacobian is singular") from error
    if not np.isfinite(result).all():
        raise ValueError("vector transformation produced a non-finite value")
    return result


def signed_permutation(matrix: object, tolerance: float = 1.0e-12) -> np.ndarray | None:
    """Return an exact cube-preserving signed permutation, or ``None``."""

    array = _finite_array(matrix, (3, 3), "spatial map")
    rounded = np.rint(array)
    if not np.allclose(array, rounded, rtol=0.0, atol=tolerance):
        return None
    integer = rounded.astype(int)
    if not np.all(np.isin(integer, (-1, 0, 1))):
        return None
    if not np.all(np.sum(np.abs(integer), axis=0) == 1):
        return None
    if not np.all(np.sum(np.abs(integer), axis=1) == 1):
        return None
    return integer


def audit_frame(frame: AffineFrame, tolerance: float = 1.0e-12) -> dict[str, object]:
    """Classify whether the current fixed cubical extractor can realize a map."""

    spatial_time_mixing = max(
        float(np.max(np.abs(frame.spatial_legs[0]))),
        float(np.max(np.abs(frame.spatial_leg_derivative[0]))),
    )
    source_time_rate = float(frame.worldline_tangent[0])
    spatial_map = frame.spatial_legs[1:, :]
    permutation = signed_permutation(spatial_map, tolerance)
    center_velocity = frame.worldline_tangent[1:]
    rotation_rate = frame.spatial_leg_derivative[1:, :]
    moving = (
        float(np.max(np.abs(center_velocity))) > tolerance
        or float(np.max(np.abs(rotation_rate))) > tolerance
    )
    ordinary_time = (
        abs(source_time_rate - 1.0) <= tolerance
        and spatial_time_mixing <= tolerance
    )
    if ordinary_time and permutation is not None and not moving:
        capability = "exact_static_cube_relabel"
        current_writer_supported = True
    elif ordinary_time and permutation is not None:
        capability = "moving_axis_aligned_ale_required"
        current_writer_supported = False
    else:
        capability = "moving_cut_surface_required"
        current_writer_supported = False
    return {
        "classification": FRAME_CONTRACT_CLASSIFICATION,
        "capability": capability,
        "current_fixed_writer_supported": current_writer_supported,
        "motional_emf_required": moving,
        "spatial_signed_permutation": (
            None if permutation is None else permutation.tolist()
        ),
        "worldline_tangent": frame.worldline_tangent.tolist(),
        "spatial_time_mixing_maximum": spatial_time_mixing,
        "spatial_leg_derivative_maximum": float(
            np.max(np.abs(frame.spatial_leg_derivative))
        ),
        "required_edge_form": (
            "pullback of F to (T,edge), reducing to (E+w_cross_B).dl when t=T"
        ),
        "required_discretization": (
            "sample moving endpoint fluxes and motional edge EMFs, then enforce the "
            "closed-surface cochain identities"
        ),
    }


class CubeSurfaceComplex:
    """Unique-edge incidence complex for an ``N x N`` cubical surface."""

    def __init__(self, cells: int):
        if not isinstance(cells, int) or cells < 1:
            raise ValueError("cube surface resolution must be a positive integer")
        self.cells = cells
        self._edge_ids: dict[tuple[int, int, int, int], int] = {}
        self.u_edge_ids: dict[str, np.ndarray] = {}
        self.u_edge_signs: dict[str, np.ndarray] = {}
        self.v_edge_ids: dict[str, np.ndarray] = {}
        self.v_edge_signs: dict[str, np.ndarray] = {}
        for name in worldtube.FACE_NAMES:
            self._build_face_edges(name)
        expected_edges = 12 * cells * cells
        if len(self._edge_ids) != expected_edges:
            raise RuntimeError(
                f"cubical surface has {len(self._edge_ids)} unique edges, "
                f"expected {expected_edges}"
            )
        self.face_edge_ids, self.face_edge_coefficients = self._build_incidence()

    @property
    def face_count(self) -> int:
        return 6 * self.cells * self.cells

    @property
    def edge_count(self) -> int:
        return len(self._edge_ids)

    def _vertex(self, name: str, u_node: int, v_node: int) -> np.ndarray:
        orientation = worldtube.ORIENTATIONS[name]
        result = np.empty(3, dtype=int)
        result[orientation.normal_axis] = (
            self.cells if orientation.normal_sign > 0 else 0
        )
        result[orientation.u_axis] = (
            u_node if orientation.u_sign > 0 else self.cells - u_node
        )
        result[orientation.v_axis] = (
            v_node if orientation.v_sign > 0 else self.cells - v_node
        )
        return result

    def _edge(self, start: np.ndarray, end: np.ndarray) -> tuple[int, int]:
        difference = end - start
        nonzero = np.flatnonzero(difference)
        if nonzero.size != 1 or abs(int(difference[nonzero[0]])) != 1:
            raise RuntimeError("cube edge endpoints are not adjacent")
        axis = int(nonzero[0])
        sign = int(difference[axis])
        lower = np.minimum(start, end)
        key = (axis, int(lower[0]), int(lower[1]), int(lower[2]))
        if key not in self._edge_ids:
            self._edge_ids[key] = len(self._edge_ids)
        return self._edge_ids[key], sign

    def _build_face_edges(self, name: str) -> None:
        cells = self.cells
        u_ids = np.empty((cells + 1, cells), dtype=np.int64)
        u_signs = np.empty_like(u_ids, dtype=np.float64)
        for v_node in range(cells + 1):
            for u in range(cells):
                edge_id, sign = self._edge(
                    self._vertex(name, u, v_node),
                    self._vertex(name, u + 1, v_node),
                )
                u_ids[v_node, u] = edge_id
                u_signs[v_node, u] = sign
        v_ids = np.empty((cells, cells + 1), dtype=np.int64)
        v_signs = np.empty_like(v_ids, dtype=np.float64)
        for v in range(cells):
            for u_node in range(cells + 1):
                edge_id, sign = self._edge(
                    self._vertex(name, u_node, v),
                    self._vertex(name, u_node, v + 1),
                )
                v_ids[v, u_node] = edge_id
                v_signs[v, u_node] = sign
        self.u_edge_ids[name] = u_ids
        self.u_edge_signs[name] = u_signs
        self.v_edge_ids[name] = v_ids
        self.v_edge_signs[name] = v_signs

    def _build_incidence(self) -> tuple[np.ndarray, np.ndarray]:
        cells = self.cells
        ids = np.empty((self.face_count, 4), dtype=np.int64)
        coefficients = np.empty((self.face_count, 4), dtype=np.float64)
        row = 0
        for name in worldtube.FACE_NAMES:
            u_ids = self.u_edge_ids[name]
            u_signs = self.u_edge_signs[name]
            v_ids = self.v_edge_ids[name]
            v_signs = self.v_edge_signs[name]
            for v in range(cells):
                for u in range(cells):
                    ids[row] = (
                        u_ids[v, u],
                        v_ids[v, u + 1],
                        u_ids[v + 1, u],
                        v_ids[v, u],
                    )
                    coefficients[row] = (
                        u_signs[v, u],
                        v_signs[v, u + 1],
                        -u_signs[v + 1, u],
                        -v_signs[v, u],
                    )
                    row += 1
        return ids, coefficients

    def curl(self, edge_values: object) -> np.ndarray:
        values = _finite_array(edge_values, (self.edge_count,), "edge cochain")
        return np.sum(
            self.face_edge_coefficients * values[self.face_edge_ids], axis=1
        )

    def transpose_curl(self, face_values: object) -> np.ndarray:
        values = _finite_array(face_values, (self.face_count,), "face cochain")
        weights = (self.face_edge_coefficients * values[:, None]).ravel()
        return np.bincount(
            self.face_edge_ids.ravel(), weights=weights, minlength=self.edge_count
        )

    def pack_edges(
        self,
        emf_u: dict[str, np.ndarray],
        emf_v: dict[str, np.ndarray],
    ) -> tuple[np.ndarray, float]:
        """Average duplicate cube-edge samples into one oriented cochain."""

        total = np.zeros(self.edge_count, dtype=np.float64)
        count = np.zeros(self.edge_count, dtype=np.int64)
        minimum = np.full(self.edge_count, np.inf)
        maximum = np.full(self.edge_count, -np.inf)
        for name in worldtube.FACE_NAMES:
            arrays = (
                (
                    _finite_array(
                        emf_u[name],
                        (self.cells + 1, self.cells),
                        f"{name} emf_u",
                    ),
                    self.u_edge_ids[name],
                    self.u_edge_signs[name],
                ),
                (
                    _finite_array(
                        emf_v[name],
                        (self.cells, self.cells + 1),
                        f"{name} emf_v",
                    ),
                    self.v_edge_ids[name],
                    self.v_edge_signs[name],
                ),
            )
            for local, edge_ids, signs in arrays:
                canonical = (signs * local).ravel()
                flat_ids = edge_ids.ravel()
                np.add.at(total, flat_ids, canonical)
                np.add.at(count, flat_ids, 1)
                np.minimum.at(minimum, flat_ids, canonical)
                np.maximum.at(maximum, flat_ids, canonical)
        if np.any(count < 1) or np.any(count > 2):
            raise RuntimeError("cubical edge occurrence count is not one or two")
        shared = count == 2
        shared_residual = (
            0.0
            if not np.any(shared)
            else float(np.max(maximum[shared] - minimum[shared]))
        )
        return total / count, shared_residual

    def unpack_edges(
        self, edge_values: object
    ) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        values = _finite_array(edge_values, (self.edge_count,), "edge cochain")
        emf_u = {
            name: self.u_edge_signs[name] * values[self.u_edge_ids[name]]
            for name in worldtube.FACE_NAMES
        }
        emf_v = {
            name: self.v_edge_signs[name] * values[self.v_edge_ids[name]]
            for name in worldtube.FACE_NAMES
        }
        return emf_u, emf_v

    def flatten_faces(self, values: dict[str, np.ndarray], name: str) -> np.ndarray:
        checked = [
            _finite_array(
                values[face], (self.cells, self.cells), f"{face} {name}"
            ).ravel()
            for face in worldtube.FACE_NAMES
        ]
        return np.concatenate(checked)


def _conjugate_gradient(
    operator: Callable[[np.ndarray], np.ndarray],
    right_hand_side: np.ndarray,
    tolerance: float,
    maximum_iterations: int,
) -> tuple[np.ndarray, int, float]:
    """Solve a zero-mean symmetric semidefinite system by projected CG."""

    rhs = np.asarray(right_hand_side, dtype=np.float64).copy()
    rhs -= np.mean(rhs)
    solution = np.zeros_like(rhs)
    residual = rhs.copy()
    direction = residual.copy()
    residual2 = float(np.dot(residual, residual))
    initial_norm = math.sqrt(residual2)
    if initial_norm == 0.0:
        return solution, 0, 0.0
    target = tolerance * initial_norm
    for iteration in range(1, maximum_iterations + 1):
        product = operator(direction)
        product -= np.mean(product)
        denominator = float(np.dot(direction, product))
        if not (math.isfinite(denominator) and denominator > 0.0):
            raise RuntimeError("surface cochain projection lost positive definiteness")
        step = residual2 / denominator
        solution += step * direction
        residual -= step * product
        residual -= np.mean(residual)
        next_residual2 = float(np.dot(residual, residual))
        norm = math.sqrt(next_residual2)
        if norm <= target:
            solution -= np.mean(solution)
            return solution, iteration, norm
        direction = residual + (next_residual2 / residual2) * direction
        residual2 = next_residual2
    raise RuntimeError(
        f"surface cochain projection did not converge after {maximum_iterations} "
        f"iterations; residual={math.sqrt(residual2):.6g}, target={target:.6g}"
    )


def project_interval_edges(
    complex_: CubeSurfaceComplex,
    target_circulation: object,
    sampled_edges: object,
    solver_tolerance: float = 2.0e-13,
    maximum_iterations: int | None = None,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Minimally correct edge samples to give an exact requested surface curl."""

    target = _finite_array(
        target_circulation, (complex_.face_count,), "target circulation"
    )
    sampled = _finite_array(sampled_edges, (complex_.edge_count,), "sampled edges")
    if not math.isfinite(solver_tolerance) or solver_tolerance <= 0.0:
        raise ValueError("solver_tolerance must be finite and positive")
    compatibility_scale = max(float(np.sum(np.abs(target))), 1.0)
    compatibility = abs(float(np.sum(target)))
    compatibility_limit = 128.0 * np.finfo(float).eps * compatibility_scale
    if compatibility > compatibility_limit:
        raise ValueError(
            "endpoint closed-surface flux changes between samples; no edge EMF "
            f"can satisfy Faraday (circulation sum={float(np.sum(target)):.6g})"
        )
    target = target - np.mean(target)
    raw_residual = target - complex_.curl(sampled)
    raw_residual -= np.mean(raw_residual)
    iterations_limit = (
        max(64, 16 * complex_.cells)
        if maximum_iterations is None
        else maximum_iterations
    )
    if not isinstance(iterations_limit, int) or iterations_limit < 1:
        raise ValueError("maximum_iterations must be a positive integer")
    potential, iterations, solver_residual = _conjugate_gradient(
        lambda values: complex_.curl(complex_.transpose_curl(values)),
        raw_residual,
        solver_tolerance,
        iterations_limit,
    )
    correction = complex_.transpose_curl(potential)
    corrected = sampled + correction
    final_residual = target - complex_.curl(corrected)
    edge_scale = max(float(np.max(np.abs(sampled))), 1.0)
    return corrected, {
        "iterations": iterations,
        "solver_l2_residual": solver_residual,
        "raw_maximum_faraday_residual": float(np.max(np.abs(raw_residual))),
        "final_maximum_faraday_residual": float(np.max(np.abs(final_residual))),
        "maximum_edge_correction": float(np.max(np.abs(correction))),
        "rms_edge_correction": float(np.sqrt(np.mean(correction**2))),
        "relative_maximum_edge_correction": float(
            np.max(np.abs(correction)) / edge_scale
        ),
    }


def project_moving_samples(
    times: object,
    faces: dict[str, worldtube.FaceData],
    solver_tolerance: float = 2.0e-13,
    maximum_iterations: int | None = None,
) -> tuple[dict[str, worldtube.FaceData], dict[str, object]]:
    """Project sampled moving-cube EMFs onto the exact CT surface complex.

    Endpoint normal fluxes and face states are retained exactly.  Edge samples
    should already contain the pulled-back motional electric field.
    """

    times_array = worldtube.validate_times(times)
    if set(faces) != set(worldtube.FACE_NAMES):
        raise ValueError("moving samples must contain exactly the six cube faces")
    checked = {
        name: worldtube.validate_face(faces[name], times_array, name)
        for name in worldtube.FACE_NAMES
    }
    reference = checked[worldtube.FACE_NAMES[0]].normal_flux.shape
    _, cells_v, cells_u = reference
    if cells_u != cells_v:
        raise ValueError("moving worldtube projection requires square face grids")
    for name in worldtube.FACE_NAMES:
        if checked[name].normal_flux.shape != reference:
            raise ValueError("moving worldtube faces must share one resolution")

    complex_ = CubeSurfaceComplex(cells_u)
    corrected_u = {
        name: np.empty_like(checked[name].emf_u) for name in worldtube.FACE_NAMES
    }
    corrected_v = {
        name: np.empty_like(checked[name].emf_v) for name in worldtube.FACE_NAMES
    }
    interval_diagnostics = []
    for interval, dt in enumerate(np.diff(times_array)):
        flux_left = {
            name: checked[name].normal_flux[interval]
            for name in worldtube.FACE_NAMES
        }
        flux_right = {
            name: checked[name].normal_flux[interval + 1]
            for name in worldtube.FACE_NAMES
        }
        target = (
            complex_.flatten_faces(flux_left, "left normal flux")
            - complex_.flatten_faces(flux_right, "right normal flux")
        ) / float(dt)
        sampled, shared_residual = complex_.pack_edges(
            {
                name: checked[name].emf_u[interval]
                for name in worldtube.FACE_NAMES
            },
            {
                name: checked[name].emf_v[interval]
                for name in worldtube.FACE_NAMES
            },
        )
        corrected, diagnostics = project_interval_edges(
            complex_,
            target,
            sampled,
            solver_tolerance=solver_tolerance,
            maximum_iterations=maximum_iterations,
        )
        unpacked_u, unpacked_v = complex_.unpack_edges(corrected)
        for name in worldtube.FACE_NAMES:
            corrected_u[name][interval] = unpacked_u[name]
            corrected_v[name][interval] = unpacked_v[name]
        diagnostics["sampled_shared_edge_residual"] = shared_residual
        diagnostics["interval"] = interval
        diagnostics["time_range"] = [
            float(times_array[interval]),
            float(times_array[interval + 1]),
        ]
        interval_diagnostics.append(diagnostics)

    result = {
        name: worldtube.FaceData(
            cell_state=checked[name].cell_state.copy(),
            normal_flux=checked[name].normal_flux.copy(),
            emf_u=corrected_u[name],
            emf_v=corrected_v[name],
        )
        for name in worldtube.FACE_NAMES
    }
    validation = worldtube.validate_worldtube(times_array, result)
    diagnostics = {
        "classification": "athenak-emri-moving-worldtube-projection-v1",
        "surface_face_cells": complex_.face_count,
        "surface_unique_edges": complex_.edge_count,
        "intervals": interval_diagnostics,
        "validation": validation,
        "warning": (
            "projection enforces topology but does not make an under-resolved moving "
            "surface physically accurate; edge corrections must converge to zero"
        ),
    }
    return result, diagnostics


def project_closed_surface_fluxes(
    times: object,
    faces: dict[str, worldtube.FaceData],
) -> tuple[dict[str, worldtube.FaceData], dict[str, object]]:
    """Minimally remove sampled magnetic monopoles at every endpoint.

    A pointwise interpolation of cell-centered magnetic data does not commute
    with the divergence operator used by the source CT mesh.  Its cut-surface
    quadrature can therefore have a small nonzero closed flux even when the
    source face field is exactly divergence-free.  On an equal-area cubical
    surface, subtracting the mean face-cell flux is the Euclidean minimum-norm
    correction satisfying the one closed-surface compatibility constraint.

    This projection is deliberately separate from ``project_moving_samples``:
    endpoint face fluxes are corrected first, and the unique edge cochain is
    then corrected to satisfy the resulting interval Faraday curls.
    """

    times_array = worldtube.validate_times(times)
    if set(faces) != set(worldtube.FACE_NAMES):
        raise ValueError("flux samples must contain exactly the six cube faces")
    checked = {
        name: worldtube.validate_face(faces[name], times_array, name)
        for name in worldtube.FACE_NAMES
    }
    reference = checked[worldtube.FACE_NAMES[0]].normal_flux.shape
    for name in worldtube.FACE_NAMES:
        if checked[name].normal_flux.shape != reference:
            raise ValueError("closed-flux projection requires one face resolution")

    corrected_flux = {
        name: checked[name].normal_flux.copy() for name in worldtube.FACE_NAMES
    }
    endpoint_diagnostics = []
    face_cell_count = sum(
        checked[name].normal_flux.shape[1]
        * checked[name].normal_flux.shape[2]
        for name in worldtube.FACE_NAMES
    )
    for endpoint in range(times_array.size):
        net_flux = float(
            sum(
                np.sum(checked[name].normal_flux[endpoint], dtype=np.float64)
                for name in worldtube.FACE_NAMES
            )
        )
        correction = -net_flux / face_cell_count
        for name in worldtube.FACE_NAMES:
            corrected_flux[name][endpoint] += correction
        final_flux = float(
            sum(
                np.sum(corrected_flux[name][endpoint], dtype=np.float64)
                for name in worldtube.FACE_NAMES
            )
        )
        flux_scale = max(
            max(
                float(np.max(np.abs(checked[name].normal_flux[endpoint])))
                for name in worldtube.FACE_NAMES
            ),
            np.finfo(float).tiny,
        )
        endpoint_diagnostics.append(
            {
                "endpoint": endpoint,
                "time": float(times_array[endpoint]),
                "raw_closed_surface_flux": net_flux,
                "final_closed_surface_flux": final_flux,
                "face_cell_correction": correction,
                "relative_maximum_flux_correction": abs(correction) / flux_scale,
            }
        )

    result = {
        name: worldtube.FaceData(
            cell_state=checked[name].cell_state.copy(),
            normal_flux=corrected_flux[name],
            emf_u=checked[name].emf_u.copy(),
            emf_v=checked[name].emf_v.copy(),
        )
        for name in worldtube.FACE_NAMES
    }
    return result, {
        "classification": "athenak-emri-closed-flux-projection-v1",
        "method": "equal-area Euclidean minimum-norm mean subtraction",
        "endpoints": endpoint_diagnostics,
        "warning": (
            "topological projection is not a substitute for spatial convergence; "
            "the relative correction must decrease as source and quadrature "
            "resolution increase"
        ),
    }


def _load_frame_contract(path: Path) -> AffineFrame:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("classification") != FRAME_CONTRACT_CLASSIFICATION:
        raise ValueError("frame contract classification is missing or unsupported")
    return AffineFrame(
        worldline_tangent=document["worldline_tangent"],
        spatial_legs=document["spatial_legs"],
        spatial_leg_derivative=document.get(
            "spatial_leg_derivative", np.zeros((4, 3))
        ),
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit_parser = subparsers.add_parser("audit")
    audit_parser.add_argument("contract", type=Path)
    project_parser = subparsers.add_parser("project-moving")
    project_parser.add_argument("input", type=Path)
    project_parser.add_argument("output", type=Path)
    project_parser.add_argument("--solver-tolerance", type=float, default=2.0e-13)
    project_parser.add_argument("--maximum-iterations", type=int)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    if arguments.command == "audit":
        diagnostics = audit_frame(_load_frame_contract(arguments.contract))
        print(json.dumps(diagnostics, indent=2))
        return
    times, faces, metadata = worldtube.read_worldtube(arguments.input)
    projected, diagnostics = project_moving_samples(
        times,
        faces,
        solver_tolerance=arguments.solver_tolerance,
        maximum_iterations=arguments.maximum_iterations,
    )
    output_metadata = dict(metadata)
    output_metadata["moving_surface_projection"] = diagnostics
    worldtube.write_worldtube(arguments.output, times, projected, output_metadata)
    print(json.dumps(diagnostics, indent=2))


if __name__ == "__main__":
    main()
