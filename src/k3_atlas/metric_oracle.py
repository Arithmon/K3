#!/usr/bin/env python3
"""Single geometric oracle for the K3 D2 Kahler ansatz.

Ambient Hessians are stored with indices ``(Z_i, conjugate(Z_j))``. Their
holomorphic pullback is therefore

    G_ab = W[:, a].T @ H @ conjugate(W[:, b]),

and scalar differentiation is ``dq = W.T @ partial_Z(q)``. No legacy
convention switch is exposed by the production API; the U-dagger contraction
exists only as an explicitly named adverse-test control.
"""
from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations_with_replacement
from typing import Iterable, Sequence

import numpy as np

METRIC_CONVENTION = "holomorphic_pullback_VT"
N_RHO = 10
MU = np.asarray([1.0, 2.0, 3.0, 5.0, 7.0, 11.0])
LAMBDA = np.vstack([np.ones(6), MU, MU ** 2])
COORD_CHARS = np.asarray([
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [0, 0, 1],
], dtype=int)


def _batch_points(z: np.ndarray) -> tuple[np.ndarray, bool]:
    array = np.asarray(z, dtype=np.complex128)
    single = array.ndim == 1
    if single:
        array = array[None, :]
    if array.ndim != 2 or array.shape[1] != 6:
        raise ValueError("Z must have shape (6,) or (n, 6)")
    return array, single


def _batch_frames(w: np.ndarray, n: int) -> tuple[np.ndarray, bool]:
    array = np.asarray(w, dtype=np.complex128)
    single = array.ndim == 2
    if single:
        array = array[None, :, :]
    if array.ndim != 3 or array.shape != (n, 6, 2):
        raise ValueError("W must have shape (6, 2) or (n, 6, 2)")
    return array, single


def q_values(z: np.ndarray) -> np.ndarray:
    points, single = _batch_points(z)
    values = np.einsum("ri,ni->nr", LAMBDA, points ** 2)
    return values[0] if single else values


def q_jacobian(z: np.ndarray) -> np.ndarray:
    points, single = _batch_points(z)
    jacobian = 2.0 * LAMBDA[None, :, :] * points[:, None, :]
    return jacobian[0] if single else jacobian


def tangent_svd(z: np.ndarray) -> np.ndarray:
    """Orthonormal horizontal holomorphic tangent frame on the unit sphere."""
    points, single = _batch_points(z)
    frames = np.empty((len(points), 6, 2), dtype=np.complex128)
    for index, point in enumerate(points):
        full = np.vstack([q_jacobian(point), point.conj()[None, :]])
        _, _, vh = np.linalg.svd(full, full_matrices=True)
        frames[index] = vh[-2:, :].conj().T
    return frames[0] if single else frames


def horizontalize_projective(z: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Normalize a projective representative and horizontalize its frame."""
    points, single = _batch_points(z)
    frames, _ = _batch_frames(w, len(points))
    scale2 = np.sum(np.abs(points) ** 2, axis=1)
    scale = np.sqrt(scale2)
    unit = points / scale[:, None]
    vertical = np.einsum("ni,nia->na", unit.conj(), frames)
    horizontal = (frames - np.einsum("ni,na->nia", unit, vertical)) \
        / scale[:, None, None]
    if single:
        return unit[0], horizontal[0]
    return unit, horizontal


def orthonormalize_holomorphic_frame(w: np.ndarray) -> np.ndarray:
    """Orthonormalize frames for the ``W.T @ H @ conjugate(W)`` convention."""
    array = np.asarray(w, dtype=np.complex128)
    single = array.ndim == 2
    if single:
        array = array[None, :, :]
    if array.ndim != 3 or array.shape[1:] != (6, 2):
        raise ValueError("W must have shape (6, 2) or (n, 6, 2)")
    gram = np.einsum("nia,nib->nab", array, array.conj())
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    if np.any(eigenvalues <= 0.0):
        raise ValueError("frame Gram matrix is not positive definite")
    transform = eigenvectors.conj() / np.sqrt(eigenvalues)[:, None, :]
    result = np.einsum("nia,nab->nib", array, transform)
    return result[0] if single else result


def radical_chart_lift(
        pivot: Sequence[int], gauge: int, u: complex, v: complex,
        signs: Sequence[int] = (1, 1, 1)) -> tuple[np.ndarray, np.ndarray]:
    """Float evaluation of a local radical chart and holomorphic frame."""
    pivot = tuple(int(i) for i in pivot)
    if len(pivot) != 3 or len(set(pivot)) != 3:
        raise ValueError("pivot must contain three distinct coordinates")
    complement = tuple(i for i in range(6) if i not in pivot)
    if gauge not in complement:
        raise ValueError("gauge coordinate must be outside the pivot")
    others = [i for i in complement if i != gauge]
    left = LAMBDA[:, list(pivot)]
    right = LAMBDA[:, list(complement)]
    reduction = -np.linalg.solve(left, right)
    order = [complement.index(gauge), complement.index(others[0]),
             complement.index(others[1])]
    reduction = reduction[:, order]
    radicands = reduction[:, 0] + reduction[:, 1] * u ** 2 \
        + reduction[:, 2] * v ** 2
    roots = np.asarray(signs, dtype=float) * np.sqrt(radicands + 0j)
    if np.min(np.abs(roots)) < 1e-10:
        raise ValueError("radical chart point is too close to a branch locus")
    point = np.zeros(6, dtype=np.complex128)
    point[gauge] = 1.0
    point[others[0]] = u
    point[others[1]] = v
    point[list(pivot)] = roots
    frame = np.zeros((6, 2), dtype=np.complex128)
    frame[others[0], 0] = 1.0
    frame[others[1], 1] = 1.0
    frame[list(pivot), 0] = reduction[:, 1] * u / roots
    frame[list(pivot), 1] = reduction[:, 2] * v / roots
    return point, frame


def build_M_with_derivatives(
        rho10: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Return ``M=L L^*`` and its ten exact parameter derivatives."""
    values = np.asarray(rho10, dtype=float)
    if values.shape != (N_RHO,):
        raise ValueError(f"rho parameters must have shape ({N_RHO},)")
    u1, ut, ua1, bar, bai, ua2, ub1, bbr, bbi, ub2 = values
    factor = np.zeros((6, 6), dtype=np.complex128)
    derivative = np.zeros((N_RHO, 6, 6), dtype=np.complex128)
    factor[0, 0] = np.exp(u1)
    derivative[0, 0, 0] = factor[0, 0]
    factor[1, 1] = np.exp(ut)
    derivative[1, 1, 1] = factor[1, 1]
    factor[2, 2] = np.exp(ua1)
    derivative[2, 2, 2] = factor[2, 2]
    factor[3, 2] = bar + 1j * bai
    derivative[3, 3, 2] = 1.0
    derivative[4, 3, 2] = 1j
    factor[3, 3] = np.exp(ua2)
    derivative[5, 3, 3] = factor[3, 3]
    factor[4, 4] = np.exp(ub1)
    derivative[6, 4, 4] = factor[4, 4]
    factor[5, 4] = bbr + 1j * bbi
    derivative[7, 5, 4] = 1.0
    derivative[8, 5, 4] = 1j
    factor[5, 5] = np.exp(ub2)
    derivative[9, 5, 5] = factor[5, 5]
    matrix = factor @ factor.conj().T
    matrix_derivative = (
        derivative @ factor.conj().T
        + factor[None, :, :] @ np.swapaxes(derivative.conj(), 1, 2))
    return matrix, matrix_derivative


def build_M(rho10: Sequence[float]) -> np.ndarray:
    return build_M_with_derivatives(rho10)[0]


def _character(indices: Iterable[int]) -> tuple[int, int, int]:
    result = np.zeros(3, dtype=int)
    for index in indices:
        result += COORD_CHARS[index]
    return tuple((result % 2).tolist())


@lru_cache(maxsize=None)
def enumerate_metric_basis(degree: int) -> tuple[dict, ...]:
    holomorphic = list(combinations_with_replacement(range(6), degree))
    result = []
    for left in holomorphic:
        for right in holomorphic:
            if _character((*left, *right)) != (0, 0, 0):
                continue
            if left < right:
                result.append({"type": "real_pair", "ij": left, "kl": right})
                result.append({"type": "imag_pair", "ij": left, "kl": right})
            elif left == right:
                result.append({"type": "self", "ij": left, "kl": right})
    return tuple(result)


@lru_cache(maxsize=1)
def metric_basis() -> tuple[dict, ...]:
    result = enumerate_metric_basis(2) + enumerate_metric_basis(3)
    if len(result) != 657:
        raise AssertionError(f"unexpected metric basis dimension {len(result)}")
    return result


def monomial_value_gradient(
        multi_index: Sequence[int], z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    points, _ = _batch_points(z)
    counts = Counter(multi_index)
    value = np.ones(len(points), dtype=np.complex128)
    for coordinate, multiplicity in counts.items():
        value *= points[:, coordinate] ** multiplicity
    gradient = np.zeros((len(points), 6), dtype=np.complex128)
    for coordinate, multiplicity in counts.items():
        derivative = np.full(len(points), complex(multiplicity))
        for other, other_multiplicity in counts.items():
            exponent = other_multiplicity - int(other == coordinate)
            if exponent:
                derivative *= points[:, other] ** exponent
        gradient[:, coordinate] = derivative
    return value, gradient


def _complex_numerator(
        z: np.ndarray, left: Sequence[int], right: Sequence[int]
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    left_value, left_gradient = monomial_value_gradient(left, z)
    right_value, right_gradient = monomial_value_gradient(right, z)
    value = left_value * right_value.conj()
    hol_gradient = left_gradient * right_value.conj()[:, None]
    raw_hessian = np.einsum(
        "ni,nj->nij", left_gradient, right_gradient.conj())
    return value, hol_gradient, raw_hessian


def real_numerator(
        z: np.ndarray, element: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    left, right = element["ij"], element["kl"]
    value, gradient, hessian = _complex_numerator(z, left, right)
    kind = element["type"]
    if kind == "self":
        return value.real, gradient, hessian
    swapped = _complex_numerator(z, right, left)
    if kind == "real_pair":
        return value.real * 2.0, gradient + swapped[1], hessian + swapped[2]
    if kind == "imag_pair":
        return -value.imag * 2.0, 1j * (gradient - swapped[1]), \
            1j * (hessian - swapped[2])
    raise ValueError(f"unknown basis element type {kind!r}")


def potential_value(z: np.ndarray, params: Sequence[float]) -> np.ndarray | float:
    points, single = _batch_points(z)
    parameters = np.asarray(params, dtype=float)
    if parameters.shape != (667,):
        raise ValueError("metric parameters must have shape (667,)")
    matrix = build_M(parameters[:N_RHO])
    rho = np.einsum("ni,ij,nj->n", points, matrix, points.conj()).real
    if np.any(rho <= 0.0):
        raise ValueError("rho is non-positive")
    scale2 = np.sum(np.abs(points) ** 2, axis=1)
    value = np.log(rho)
    for coefficient, element in zip(parameters[N_RHO:], metric_basis()):
        if coefficient == 0.0:
            continue
        numerator, _, _ = real_numerator(points, element)
        degree = len(element["ij"])
        value += coefficient * numerator * scale2 ** (-degree)
    return float(value[0]) if single else value


def fs_ambient_hessian(z: np.ndarray) -> np.ndarray:
    points, single = _batch_points(z)
    scale2 = np.sum(np.abs(points) ** 2, axis=1)
    identity = np.eye(6, dtype=np.complex128)[None, :, :]
    hessian = identity / scale2[:, None, None] - np.einsum(
        "ni,nj->nij", points.conj(), points) / scale2[:, None, None] ** 2
    return hessian[0] if single else hessian


def ambient_hessian_full(z: np.ndarray, params: Sequence[float]) -> np.ndarray:
    """Full ambient mixed Hessian, valid away from any chosen gauge/frame."""
    points, single = _batch_points(z)
    parameters = np.asarray(params, dtype=float)
    if parameters.shape != (667,):
        raise ValueError("metric parameters must have shape (667,)")
    matrix = build_M(parameters[:N_RHO])
    matrix_zbar = points.conj() @ matrix.T
    rho = np.einsum("ni,ni->n", points, matrix_zbar).real
    hessian = matrix[None, :, :] / rho[:, None, None] - np.einsum(
        "ni,nj->nij", matrix_zbar, matrix_zbar.conj()) \
        / rho[:, None, None] ** 2
    scale2 = np.sum(np.abs(points) ** 2, axis=1)
    identity = np.eye(6, dtype=np.complex128)[None, :, :]
    for coefficient, element in zip(parameters[N_RHO:], metric_basis()):
        if coefficient == 0.0:
            continue
        numerator, hol_gradient, raw_hessian = real_numerator(points, element)
        degree = len(element["ij"])
        anti_gradient = hol_gradient.conj()
        sd = scale2 ** (-degree)
        term = raw_hessian * sd[:, None, None]
        term -= degree * np.einsum(
            "ni,nj->nij", hol_gradient, points) \
            * scale2[:, None, None] ** (-degree - 1)
        term -= degree * np.einsum(
            "ni,nj->nij", points.conj(), anti_gradient) \
            * scale2[:, None, None] ** (-degree - 1)
        term += degree * (degree + 1) * numerator[:, None, None] \
            * np.einsum("ni,nj->nij", points.conj(), points) \
            * scale2[:, None, None] ** (-degree - 2)
        term -= degree * numerator[:, None, None] * identity \
            * scale2[:, None, None] ** (-degree - 1)
        hessian += coefficient * term
    return hessian[0] if single else hessian


def pullback_holomorphic(hessian: np.ndarray, w: np.ndarray) -> np.ndarray:
    h = np.asarray(hessian, dtype=np.complex128)
    single = h.ndim == 2
    if single:
        h = h[None, :, :]
    frames, _ = _batch_frames(w, len(h))
    metric = np.einsum("nia,nij,njb->nab", frames, h, frames.conj())
    return metric[0] if single else metric


def pullback_legacy_udagger(hessian: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Retracted control used only to prove that adverse tests can fail."""
    h = np.asarray(hessian, dtype=np.complex128)
    single = h.ndim == 2
    if single:
        h = h[None, :, :]
    frames, _ = _batch_frames(w, len(h))
    metric = np.einsum("nia,nij,njb->nab", frames.conj(), h, frames)
    return metric[0] if single else metric


def metric_full(z: np.ndarray, w: np.ndarray, params: Sequence[float]) -> np.ndarray:
    return pullback_holomorphic(ambient_hessian_full(z, params), w)


def metric_full_legacy_control(
        z: np.ndarray, w: np.ndarray, params: Sequence[float]) -> np.ndarray:
    return pullback_legacy_udagger(ambient_hessian_full(z, params), w)


def fs_metric(z: np.ndarray, w: np.ndarray) -> np.ndarray:
    return pullback_holomorphic(fs_ambient_hessian(z), w)


def feature_values_and_dq(
        z: np.ndarray, w: np.ndarray, basis: Sequence[dict], *,
        legacy_control: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate real projective features and their true chain derivatives."""
    points, _ = _batch_points(z)
    frames, _ = _batch_frames(w, len(points))
    scale2 = np.sum(np.abs(points) ** 2, axis=1)
    ds = np.einsum("ni,nia->na", points.conj(), frames)
    values = np.empty((len(points), len(basis)), dtype=float)
    derivatives = np.empty((len(points), len(basis), 2), dtype=np.complex128)
    for index, element in enumerate(basis):
        numerator, hol_gradient, _ = real_numerator(points, element)
        degree = len(element["ij"])
        if legacy_control:
            dnumerator = np.einsum(
                "nia,ni->na", frames.conj(), hol_gradient)
        else:
            dnumerator = np.einsum("nia,ni->na", frames, hol_gradient)
        sd = scale2 ** (-degree)
        values[:, index] = numerator * sd
        derivatives[:, index] = sd[:, None] * (
            dnumerator - degree * numerator[:, None] * ds / scale2[:, None])
    return values, derivatives


def metric_basis_cache(
        z: np.ndarray, v: np.ndarray,
        basis: Sequence[dict] | None = None) -> dict[str, np.ndarray]:
    """Correct coefficient-linear cache used by the v2 optimizer."""
    points, _ = _batch_points(z)
    frames, _ = _batch_frames(v, len(points))
    elements = tuple(metric_basis() if basis is None else basis)
    scale2 = np.sum(np.abs(points) ** 2, axis=1)
    ds = np.einsum("ni,nia->na", points.conj(), frames)
    gram = np.einsum("nia,nib->nab", frames, frames.conj())
    tensors = np.empty((len(points), len(elements), 2, 2), dtype=np.complex128)
    values = np.empty((len(points), len(elements)), dtype=float)
    for index, element in enumerate(elements):
        numerator, hol_gradient, raw_hessian = real_numerator(points, element)
        degree = len(element["ij"])
        projected_raw = pullback_holomorphic(raw_hessian, frames)
        dnumerator = np.einsum("nia,ni->na", frames, hol_gradient)
        sd = scale2 ** (-degree)
        tensor = projected_raw * sd[:, None, None]
        tensor -= degree * np.einsum(
            "na,nb->nab", dnumerator, ds.conj()) \
            * scale2[:, None, None] ** (-degree - 1)
        tensor -= degree * np.einsum(
            "na,nb->nab", ds, dnumerator.conj()) \
            * scale2[:, None, None] ** (-degree - 1)
        tensor += degree * (degree + 1) * numerator[:, None, None] \
            * np.einsum("na,nb->nab", ds, ds.conj()) \
            * scale2[:, None, None] ** (-degree - 2)
        tensor -= degree * numerator[:, None, None] * gram \
            * scale2[:, None, None] ** (-degree - 1)
        tensors[:, index] = tensor
        values[:, index] = numerator * sd
    return {
        "tensors": tensors,
        "values": values,
        "degrees": np.asarray([len(element["ij"]) for element in elements]),
    }


def rho_metric(z: np.ndarray, v: np.ndarray, rho10: Sequence[float]) -> np.ndarray:
    points, single = _batch_points(z)
    frames, _ = _batch_frames(v, len(points))
    matrix = build_M(rho10)
    matrix_zbar = points.conj() @ matrix.T
    rho = np.einsum("ni,ni->n", points, matrix_zbar).real
    first = np.einsum("nia,ij,njb->nab", frames, matrix, frames.conj())
    directional = np.einsum("nia,ni->na", frames, matrix_zbar)
    result = first / rho[:, None, None] - np.einsum(
        "na,nb->nab", directional, directional.conj()) \
        / rho[:, None, None] ** 2
    return result[0] if single else result


def rho_metric_with_derivatives(
        z: np.ndarray, v: np.ndarray,
        rho10: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    """Return the rho metric and exact derivatives in the ten rho slots."""
    points, single = _batch_points(z)
    frames, _ = _batch_frames(v, len(points))
    matrix, matrix_derivative = build_M_with_derivatives(rho10)
    matrix_zbar = points.conj() @ matrix.T
    matrix_zbar_derivative = np.einsum(
        "kij,nj->nki", matrix_derivative, points.conj())
    rho = np.einsum("ni,ni->n", points, matrix_zbar).real
    rho_derivative = np.einsum(
        "ni,nki->nk", points, matrix_zbar_derivative).real
    first = np.einsum("nia,ij,njb->nab", frames, matrix, frames.conj())
    first_derivative = np.einsum(
        "nia,kij,njb->nkab", frames, matrix_derivative, frames.conj())
    directional = np.einsum("nia,ni->na", frames, matrix_zbar)
    directional_derivative = np.einsum(
        "nia,nki->nka", frames, matrix_zbar_derivative)
    outer = np.einsum("na,nb->nab", directional, directional.conj())
    outer_derivative = (
        np.einsum("nka,nb->nkab", directional_derivative,
                  directional.conj())
        + np.einsum("na,nkb->nkab", directional,
                    directional_derivative.conj()))
    metric = first / rho[:, None, None] - outer / rho[:, None, None] ** 2
    metric_derivative = (
        first_derivative / rho[:, None, None, None]
        - first[:, None, :, :] * rho_derivative[:, :, None, None]
        / rho[:, None, None, None] ** 2
        - outer_derivative / rho[:, None, None, None] ** 2
        + 2.0 * outer[:, None, :, :] * rho_derivative[:, :, None, None]
        / rho[:, None, None, None] ** 3)
    if single:
        return metric[0], metric_derivative[0]
    return metric, metric_derivative
