#!/usr/bin/env python3
"""Versioned witness loading with an explicit retraction boundary.

The active K3 metric datum is never inferred from a historical filename.
Retracted byte streams remain readable only through an explicit audit opt-in.

ACTIVE DATUM (reconciled 2026-07-16)
------------------------------------
The active witness is the frozen witness manifest `results/k3_closedform_witness_kahler_v2.npz`
(schema `k3_kahler_witness_v2`: native coeffs218 + M, gauge det M = 1), produced
by `witness_manifest.py` and certified by the witness check `witness_check.py`
(11/11 PASS). Its bytes are pinned by `artifact_sha256` in the sidecar manifest
`results/witness_manifest.json`; the frozen npz is NEVER rewritten — this
module validates it against the sidecar instead of stamping metadata into it.

Convention: the R3 datum is built on `kahler_metric.py`, whose chart
metric is G = del delbar K~ with the HOLOMORPHIC chain-rule contraction (V^T),
i.e. exactly the `holomorphic_pullback_VT` convention required by
`results/retracted/RETRACTED.json` after the 2026-07-13 retraction of the
U-dagger contraction. See notes the metric convention retraction
and the witness frozen manifest check.

LEGACY VT-667 LINE
------------------
The a review VT refit line (schema `k3_cy_witness_v2`, params_full(667)) ended in
NO-GO on 2026-07-13 (legacy_fit_refusal): no active witness ever
existed under that schema. Its loader path is kept for explicit-path audit
replays used by the historical d2_vt scripts.
"""
from __future__ import annotations

import hashlib
import io
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

ROOT = Path(__file__).resolve().parent

ACTIVE_WITNESS_V2 = ROOT / "data" / "k3_closedform_witness_kahler_v2.npz"
ACTIVE_MANIFEST_V2 = ROOT / "data" / "witness_manifest.json"
RETRACTION_REGISTRY = ROOT / "data" / "retracted" / "RETRACTED.json"
REQUIRED_METRIC_CONVENTION = "holomorphic_pullback_VT"
REQUIRED_WITNESS_SCHEMA = "k3_kahler_witness_v2"

# Legacy VT-667 line (no-go 2026-07-13, audit replays only).
LEGACY_VT_SCHEMA = "k3_cy_witness_v2"
LEGACY_VT_PARAMS_V2 = ROOT / "data" / "legacy_fit_parameters.npz"


class WitnessArtifactError(RuntimeError):
    """Raised when a witness is missing, retracted or incorrectly versioned."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def retraction_registry() -> dict[str, Any]:
    return json.loads(RETRACTION_REGISTRY.read_text(encoding="utf-8"))


def _retracted_by_sha() -> dict[str, dict[str, Any]]:
    return {entry["sha256"]: entry
            for entry in retraction_registry()["artifacts"]}


def _scalar_text(payload: dict[str, np.ndarray], key: str) -> str:
    if key not in payload:
        raise WitnessArtifactError(f"witness is missing required field {key!r}")
    value = payload[key]
    if value.shape != ():
        raise WitnessArtifactError(f"witness field {key!r} must be scalar")
    return str(value.item())


def _validate_sha_text(value: str, field: str) -> None:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise WitnessArtifactError(f"witness field {field!r} is not a SHA256")


def _validate_kahler_v2_payload(
        payload: dict[str, np.ndarray], digest: str) -> dict[str, Any]:
    """Validate the R3 kahler witness against its embedded + sidecar manifest.

    Returns the parsed manifest. The frozen npz carries no schema fields by
    design (its bytes were pinned by the R3 check before this registry knew the
    schema); integrity is enforced through the sidecar pin instead.
    """
    for key, shape in (("p9", (9,)), ("c208", (208,)), ("coeffs218", (218,))):
        array = payload.get(key)
        if array is None or array.shape != shape:
            raise WitnessArtifactError(
                f"kahler witness field {key!r} must have shape {shape}")
        if not np.isfinite(array).all():
            raise WitnessArtifactError(f"kahler witness field {key!r} is not finite")
    matrix = payload.get("M")
    if matrix is None or matrix.shape != (6, 6) or matrix.dtype.kind != "c":
        raise WitnessArtifactError("kahler witness M must be complex (6, 6)")
    if not np.isfinite(matrix).all().item():
        raise WitnessArtifactError("kahler witness M is not finite")
    determinant = complex(np.linalg.det(matrix))
    if abs(determinant - 1.0) > 1e-9:
        raise WitnessArtifactError(
            f"kahler witness gauge det M = 1 violated: det M = {determinant}")

    manifest = json.loads(_scalar_text(payload, "manifest"))
    provenance = manifest["provenance"]
    for field in ("hash_C", "hash_basis_B3", "hash_fit_vector"):
        value = _scalar_text(payload, field)
        _validate_sha_text(value, field)
        if value != provenance[field]:
            raise WitnessArtifactError(
                f"witness field {field!r} disagrees with embedded manifest")
    protocol = manifest["protocol"]
    if int(payload["seed_fit"].item()) != protocol["seed_fit"]:
        raise WitnessArtifactError("seed_fit disagrees with embedded manifest")
    if int(payload["n_base_fit"].item()) != protocol["n_base_fit"]:
        raise WitnessArtifactError("n_base_fit disagrees with embedded manifest")
    if float(payload["var_fit"].item()) != manifest["fit"]["var_fit"]:
        raise WitnessArtifactError("var_fit disagrees with embedded manifest")

    if not ACTIVE_MANIFEST_V2.exists():
        raise WitnessArtifactError(
            f"sidecar manifest does not exist: {ACTIVE_MANIFEST_V2}")
    sidecar = json.loads(ACTIVE_MANIFEST_V2.read_text(encoding="utf-8"))
    pinned = sidecar.pop("artifact_sha256", None)
    if pinned != digest:
        raise WitnessArtifactError(
            "witness bytes do not match the frozen witness manifest pin "
            f"(sidecar {pinned}, artifact {digest})")
    if sidecar != manifest:
        raise WitnessArtifactError(
            "embedded manifest disagrees with sidecar witness_manifest.json")
    return manifest


def _validate_legacy_vt_payload(payload: dict[str, np.ndarray]) -> None:
    parameters = payload.get("params_full")
    rho = payload.get("rho10")
    phi2 = payload.get("phi2")
    phi3 = payload.get("phi3")
    if parameters is None or parameters.shape != (667,):
        raise WitnessArtifactError("active witness params_full must have shape (667,)")
    if rho is None or rho.shape != (10,):
        raise WitnessArtifactError("active witness rho10 must have shape (10,)")
    if phi2 is None or phi2.shape != (97,):
        raise WitnessArtifactError("active witness phi2 must have shape (97,)")
    if phi3 is None or phi3.shape != (560,):
        raise WitnessArtifactError("active witness phi3 must have shape (560,)")
    if not np.isfinite(parameters).all():
        raise WitnessArtifactError("active witness contains non-finite parameters")
    if not np.array_equal(parameters[:10], rho):
        raise WitnessArtifactError("rho10 does not match params_full")
    if not np.array_equal(parameters[10:107], phi2):
        raise WitnessArtifactError("phi2 does not match params_full")
    if not np.array_equal(parameters[107:], phi3):
        raise WitnessArtifactError("phi3 does not match params_full")
    test_z = payload.get("test_Z")
    test_logjjh = payload.get("test_logjjh")
    if test_z is None or test_z.ndim != 2 or test_z.shape[1] != 6:
        raise WitnessArtifactError("active witness test_Z must have shape (n, 6)")
    if test_logjjh is None or test_logjjh.shape != (len(test_z),):
        raise WitnessArtifactError("test_logjjh must match test_Z")
    for field in ("source_params_sha256", "sample_parent_sha256"):
        _validate_sha_text(_scalar_text(payload, field), field)
    commit = _scalar_text(payload, "git_commit")
    if len(commit) != 40 or any(character not in "0123456789abcdef" for character in commit):
        raise WitnessArtifactError("witness git_commit must be a full hexadecimal SHA")
    _scalar_text(payload, "created_utc")
    schema = _scalar_text(payload, "schema_version")
    status = _scalar_text(payload, "artifact_status")
    convention = _scalar_text(payload, "metric_convention")
    if schema != LEGACY_VT_SCHEMA:
        raise WitnessArtifactError(
            f"unsupported witness schema {schema!r}; expected {LEGACY_VT_SCHEMA!r}")
    if status != "ACTIVE":
        raise WitnessArtifactError(
            f"refusing witness with artifact_status={status!r}")
    if convention != REQUIRED_METRIC_CONVENTION:
        raise WitnessArtifactError(
            f"metric convention {convention!r} is not "
            f"{REQUIRED_METRIC_CONVENTION!r}")


def load_witness_artifact(
        path: Path | str | None = None, *, allow_retracted: bool = False,
        require_active_metadata: bool = True) -> dict[str, np.ndarray]:
    """Load a witness while enforcing retraction and convention metadata.

    With no ``path``, loads and validates the ACTIVE R3 kahler witness
    (`coeffs218` datum). ``allow_retracted`` is intended only for adverse
    tests and historical replay. It never upgrades a retracted artifact to
    active status.
    """
    artifact = Path(path) if path is not None else ACTIVE_WITNESS_V2
    if not artifact.is_absolute():
        cwd_candidate = (Path.cwd() / artifact).resolve()
        repository_candidate = (ROOT.parent / artifact).resolve()
        artifact = cwd_candidate if cwd_candidate.exists() else repository_candidate
    if not artifact.exists():
        if artifact == ACTIVE_WITNESS_V2:
            raise WitnessArtifactError(
                "no active K3 witness v2 exists; D2 refit/frozen manifest is required")
        raise WitnessArtifactError(f"witness artifact does not exist: {artifact}")

    digest = sha256_file(artifact)
    retracted = _retracted_by_sha().get(digest)
    if retracted is not None and not allow_retracted:
        raise WitnessArtifactError(
            f"refusing retracted artifact {artifact} ({digest}); "
            "pass allow_retracted=True only for an explicit audit replay")

    with np.load(artifact, allow_pickle=False) as archive:
        payload = {key: archive[key] for key in archive.files}

    if retracted is not None:
        payload["artifact_status"] = np.asarray("RETRACTED")
        payload["artifact_sha256"] = np.asarray(digest)
        payload["archived_path"] = np.asarray(retracted["archived_path"])
        return payload

    if require_active_metadata:
        if "coeffs218" in payload:
            _validate_kahler_v2_payload(payload, digest)
            payload["schema_version"] = np.asarray(REQUIRED_WITNESS_SCHEMA)
            payload["artifact_status"] = np.asarray("ACTIVE")
            payload["metric_convention"] = np.asarray(REQUIRED_METRIC_CONVENTION)
        else:
            _validate_legacy_vt_payload(payload)
    payload["artifact_sha256"] = np.asarray(digest)
    return payload


def load_active_witness() -> dict[str, np.ndarray]:
    """Load the ACTIVE R3 kahler witness (the only entry point R4 should use)."""
    return load_witness_artifact()


def load_canonical_MH() -> dict[str, Any]:
    """C47 (review GPT P0a2-A) : expose the canonical hermitian datum.

    The serialized witness M is NOT bit-hermitian (anti-hermitian residual
    ~1.45e-23, caught by C37). The analytic object every metric engine
    consumes is H = (M+M†)/2 — canonize it HERE, once, and let consumers
    request `M_H_canonical` and check on `sha256_MH` instead of re-deriving
    it locally. `M_serialized` stays available for byte-level audits."""
    payload = load_witness_artifact()
    M = np.asarray(payload["M"], complex)
    M_H = 0.5 * (M + M.conj().T)
    if not np.array_equal(M_H, M_H.conj().T):
        raise WitnessArtifactError("canonical M_H is not bit-hermitian")
    return {
        "M_serialized": M,
        "M_H_canonical": M_H,
        "sha256_MH": hashlib.sha256(
            np.ascontiguousarray(M_H).tobytes()).hexdigest(),
        "antihermitian_residual": float(np.max(np.abs(M - M_H))),
        "witness_sha256": str(payload["artifact_sha256"].item()),
        "coeffs218": np.asarray(payload["coeffs218"], float)}


def assert_no_active_v1_path() -> None:
    historical = ROOT / "notebooks" / "k3_closedform_witness_kahler_v1.npz"
    if historical.exists():
        raise WitnessArtifactError(
            f"retracted witness unexpectedly restored at active path {historical}")


def _selftest() -> int:
    import shutil
    import tempfile

    failures: list[str] = []

    def check(name: str, condition: bool, detail: str = "") -> None:
        status = "PASS" if condition else "FAIL"
        print(f"[{status}] {name}" + (f" — {detail}" if detail else ""))
        if not condition:
            failures.append(name)

    assert_no_active_v1_path()
    check("S0 no v1 at active path", True)

    payload = load_active_witness()
    check("S1 active witness loads",
          str(payload["artifact_status"].item()) == "ACTIVE",
          str(payload["artifact_sha256"].item())[:16])
    check("S2 schema + convention",
          str(payload["schema_version"].item()) == REQUIRED_WITNESS_SCHEMA
          and str(payload["metric_convention"].item()) == REQUIRED_METRIC_CONVENTION)
    check("S3 native datum shape", payload["coeffs218"].shape == (218,))

    with tempfile.TemporaryDirectory() as tmp:
        tampered_path = Path(tmp) / "tampered.npz"
        arrays = {key: np.array(value) for key, value in payload.items()
                  if key not in ("artifact_status", "artifact_sha256",
                                 "schema_version", "metric_convention")}
        arrays["coeffs218"] = arrays["coeffs218"].copy()
        arrays["coeffs218"][0] += 1e-9
        np.savez(tampered_path, **arrays)
        try:
            load_witness_artifact(tampered_path)
        except WitnessArtifactError as exc:
            check("S4 tampered bytes refused", "pin" in str(exc), str(exc)[:70])
        else:
            check("S4 tampered bytes refused", False)

        missing_path = Path(tmp) / "missing_manifest.npz"
        stripped = {key: value for key, value in arrays.items() if key != "manifest"}
        np.savez(missing_path, **stripped)
        try:
            load_witness_artifact(missing_path)
        except WitnessArtifactError:
            check("S5 missing manifest refused", True)
        else:
            check("S5 missing manifest refused", False)

        moved_path = Path(tmp) / "moved.npz"
        shutil.copy2(ACTIVE_WITNESS_V2, moved_path)
        moved = load_witness_artifact(moved_path)
        check("S6 byte-identical copy accepted",
              str(moved["artifact_sha256"].item())
              == str(payload["artifact_sha256"].item()))

    retracted_v1 = ROOT / "data" / "retracted" / "k3_closedform_witness_kahler_v1.npz"
    try:
        load_witness_artifact(retracted_v1)
    except WitnessArtifactError as exc:
        check("S7 retracted v1 refused", "retracted" in str(exc))
    else:
        check("S7 retracted v1 refused", False)
    audit = load_witness_artifact(retracted_v1, allow_retracted=True)
    check("S8 audit opt-in stays RETRACTED",
          str(audit["artifact_status"].item()) == "RETRACTED")

    print(f"selftest: {8 - len(failures) + 1}/9 PASS"
          if not failures else f"selftest FAILURES: {failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    assert_no_active_v1_path()
    if len(sys.argv) == 2 and sys.argv[1] == "--selftest":
        raise SystemExit(_selftest())
    if len(sys.argv) == 3 and sys.argv[1] == "--audit-retracted":
        audit_path = Path(sys.argv[2])
        payload = load_witness_artifact(audit_path, allow_retracted=True)
        print(json.dumps({
            "status": str(payload["artifact_status"].item()),
            "sha256": str(payload["artifact_sha256"].item()),
            "path": str(audit_path),
            "mode": "EXPLICIT_AUDIT_ONLY",
        }, indent=2))
    elif len(sys.argv) == 1:
        try:
            payload = load_witness_artifact()
        except WitnessArtifactError as exc:
            print(f"ACTIVE_WITNESS_UNAVAILABLE: {exc}")
            raise SystemExit(1)
        manifest = json.loads(str(payload["manifest"].item()))
        print(f"ACTIVE_WITNESS: {ACTIVE_WITNESS_V2}")
        print(f"  schema:     {payload['schema_version'].item()}")
        print(f"  convention: {payload['metric_convention'].item()}")
        print(f"  sha256:     {payload['artifact_sha256'].item()}")
        print(f"  var_fit:    {manifest['fit']['var_fit']:.6g} "
              f"(seed {manifest['protocol']['seed_fit']}, "
              f"n_base {manifest['protocol']['n_base_fit']})")
    else:
        raise SystemExit(
            "usage: witness_registry.py [--selftest | --audit-retracted PATH]")
