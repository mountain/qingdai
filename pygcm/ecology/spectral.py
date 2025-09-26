from __future__ import annotations

import os
import numpy as np
from dataclasses import dataclass


@dataclass
class SpectralBands:
    """
    Discrete spectral bands for visible/near-visible shortwave.
    Wavelength unit: nm.
    """
    nbands: int
    lambda_edges: np.ndarray  # shape [NB+1]
    lambda_centers: np.ndarray  # shape [NB]
    delta_lambda: np.ndarray  # shape [NB]

    def as_tuple(self):
        return self.nbands, self.lambda_edges, self.lambda_centers, self.delta_lambda


def make_bands(nbands: int | None = None,
               lam0_nm: float | None = None,
               lam1_nm: float | None = None) -> SpectralBands:
    """
    Construct equally spaced spectral bands in [lam0, lam1] (nm).
    Defaults: NB from env QD_ECO_SPECTRAL_BANDS (fallback 16),
              range from env QD_ECO_SPECTRAL_RANGE_NM (fallback 380,780).
    """
    if nbands is None:
        try:
            nbands = int(os.getenv("QD_ECO_SPECTRAL_BANDS", "16"))
        except Exception:
            nbands = 16
    if lam0_nm is None or lam1_nm is None:
        rng = os.getenv("QD_ECO_SPECTRAL_RANGE_NM", "380,780")
        try:
            lam0_nm, lam1_nm = [float(x.strip()) for x in rng.split(",")]
        except Exception:
            lam0_nm, lam1_nm = 380.0, 780.0
    # Guard
    nbands = max(1, int(nbands))
    lam0_nm = float(lam0_nm)
    lam1_nm = float(lam1_nm)
    if lam1_nm <= lam0_nm:
        lam0_nm, lam1_nm = 380.0, 780.0

    edges = np.linspace(lam0_nm, lam1_nm, nbands + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = (edges[1:] - edges[:-1])
    return SpectralBands(nbands=nbands,
                         lambda_edges=edges.astype(float),
                         lambda_centers=centers.astype(float),
                         delta_lambda=widths.astype(float))


def _rayleigh_weight(centers_nm: np.ndarray,
                     t0: float,
                     lref_nm: float,
                     eta: float) -> np.ndarray:
    """
    Simplified Rayleigh transmittance weight ~ T0 * (λ/λ_ref)^eta (eta≈4).
    We will normalize later, so absolute scale is not critical.
    """
    lam = np.maximum(1e-6, centers_nm)
    w = t0 * (lam / max(1e-6, lref_nm)) ** float(eta)
    w = np.clip(w, 0.0, None)
    return w


def _greenish_leaf_reflectance(centers_nm: np.ndarray) -> np.ndarray:
    """
    Minimalistic green-ish leaf reflectance template in [0,1] for M1:
    - Modest baseline reflectance 0.25
    - Gaussian bump centered near 550 nm (green) with width ~60 nm and height 0.15
    Clamp to [0,1].
    """
    mu = 550.0
    sigma = 60.0
    base = 0.25
    bump = 0.15 * np.exp(-((centers_nm - mu) ** 2) / (2.0 * sigma ** 2))
    r = np.clip(base + bump, 0.0, 1.0)
    return r


def toa_to_surface_bands(I_total: np.ndarray,
                         cloud_eff: np.ndarray | float,
                         bands: SpectralBands,
                         mode: str | None = None) -> np.ndarray:
    """
    Convert total shortwave at TOA/surface proxy (I_total, W/m^2) into band-averaged intensities.
    For M1, we compute a global band weight vector and apply to each grid cell:
      I_b = (w_b / sum(w_b)) * I_total

    Args:
        I_total: 2D array (lat x lon) current step insolation proxy (gcm.isr = insA + insB)
        cloud_eff: unused in M1 weighting (kept for future extension)
        bands: SpectralBands
        mode: 'simple' or 'rayleigh' (default from env QD_ECO_TOA_TO_SURF_MODE; fallback 'simple')

    Returns:
        I_b: 3D array [NB, n_lat, n_lon] with band-averaged intensities, sum over bands ~= I_total
    """
    nlat, nlon = I_total.shape
    mode = (mode or os.getenv("QD_ECO_TOA_TO_SURF_MODE", "simple")).strip().lower()

    if mode == "rayleigh":
        try:
            t0 = float(os.getenv("QD_ECO_RAYLEIGH_T0", "0.9"))
        except Exception:
            t0 = 0.9
        try:
            lref = float(os.getenv("QD_ECO_RAYLEIGH_LREF_NM", "550"))
        except Exception:
            lref = 550.0
        try:
            eta = float(os.getenv("QD_ECO_RAYLEIGH_ETA", "4.0"))
        except Exception:
            eta = 4.0
        w = _rayleigh_weight(bands.lambda_centers, t0, lref, eta)
    else:
        # simple: flat weighting across bands
        w = np.ones_like(bands.lambda_centers, dtype=float)

    # Normalize to unit sum
    wsum = float(np.sum(w)) + 1e-12
    wn = w / wsum

    # Allocate and broadcast
    I_b = np.empty((bands.nbands, nlat, nlon), dtype=float)
    for k in range(bands.nbands):
        I_b[k, :, :] = wn[k] * I_total
    return I_b


def band_weights_from_mode(bands: SpectralBands,
                           mode: str | None = None) -> np.ndarray:
    """
    Return normalized band weights (sum=1) for given mode (simple|rayleigh).
    Useful for reducing band albedo to scalar alpha by alpha = sum(A_b * w_b).
    """
    mode = (mode or os.getenv("QD_ECO_TOA_TO_SURF_MODE", "simple")).strip().lower()
    if mode == "rayleigh":
        try:
            t0 = float(os.getenv("QD_ECO_RAYLEIGH_T0", "0.9"))
        except Exception:
            t0 = 0.9
        try:
            lref = float(os.getenv("QD_ECO_RAYLEIGH_LREF_NM", "550"))
        except Exception:
            lref = 550.0
        try:
            eta = float(os.getenv("QD_ECO_RAYLEIGH_ETA", "4.0"))
        except Exception:
            eta = 4.0
        w = _rayleigh_weight(bands.lambda_centers, t0, lref, eta)
    else:
        w = np.ones_like(bands.lambda_centers, dtype=float)
    wsum = float(np.sum(w)) + 1e-12
    return w / wsum


def default_leaf_reflectance(bands: SpectralBands) -> np.ndarray:
    """
    Return a default leaf reflectance curve R_leaf[NB] in [0,1] for M1.
    """
    return _greenish_leaf_reflectance(bands.lambda_centers)


def absorbance_from_genes(bands: SpectralBands, genes) -> np.ndarray:
    """
    Compute band absorbance A_b[NB] in [0,1] from a gene's absorption_peaks definition.

    Expected genes.absorption_peaks: iterable of peaks where each peak has
      - center_nm: float (peak center wavelength, nm)
      - width_nm:  float (Gaussian width parameter, nm; interpreted as sigma; if FWHM is provided upstream,
                           consider converting before; here we assume sigma directly)
      - height:    float in [0,1] (peak amplitude)
    The final absorbance is clipped to [0,1]. If no peaks provided, fall back to:
      A_b = 1 - default_leaf_reflectance(bands)
    """
    NB = int(getattr(bands, "nbands", 1))
    lam = np.asarray(getattr(bands, "lambda_centers", np.linspace(400.0, 700.0, NB)), dtype=float).ravel()
    if lam.shape[0] != NB:
        lam = np.linspace(float(lam.min(initial=400.0)), float(lam.max(initial=700.0)), NB)

    # Try read peaks from genes
    peaks = []
    try:
        peaks = getattr(genes, "absorption_peaks", []) or []
    except Exception:
        peaks = []

    if not peaks:
        # Fallback: absorption is 1 - default reflectance
        R_leaf = default_leaf_reflectance(bands)
        return np.clip(1.0 - np.asarray(R_leaf, dtype=float).ravel(), 0.0, 1.0)

    A = np.zeros((NB,), dtype=float)

    def _get_attr(obj, name: str, default: float) -> float:
        # support dataclass-like, SimpleNamespace, dict, or tuple/list
        if hasattr(obj, name):
            try:
                return float(getattr(obj, name))
            except Exception:
                pass
        if isinstance(obj, dict) and name in obj:
            try:
                return float(obj[name])
            except Exception:
                pass
        return float(default)

    for pk in peaks:
        c = _get_attr(pk, "center_nm", 550.0)
        w = _get_attr(pk, "width_nm", 50.0)
        h = _get_attr(pk, "height", 0.5)
        # Guard values
        w = max(1e-3, float(w))
        h = float(np.clip(h, 0.0, 1.0))
        # Gaussian peak (sigma = width_nm)
        A += h * np.exp(-((lam - c) ** 2) / (2.0 * w ** 2))

    # Clip to [0,1]
    A = np.clip(A, 0.0, 1.0)
    return A


# --------- Spectral physics: main-sequence effective temperatures and blackbody bands ---------

_T_SUN = 5778.0  # K
_h = 6.62607015e-34  # J*s
_c = 2.99792458e8    # m/s
_kB = 1.380649e-23   # J/K

def estimate_teff_from_LM(L_ratio: float, M_ratio: float, j: float = 0.8, T_sun: float = _T_SUN) -> float:
    """
    Estimate main-sequence effective temperature (K) from luminosity and mass ratios using:
        T = T_sun * (L/L_sun)^(1/4) * (M/M_sun)^(-j/2)
    j in [0.5..1.0] depending on mass regime; defaults to 0.8 for 0.5–2 Msun.
    """
    L_ratio = float(max(L_ratio, 1e-12))
    M_ratio = float(max(M_ratio, 1e-12))
    return float(T_sun * (L_ratio ** 0.25) * (M_ratio ** (-0.5 * j)))


def _planck_lambda_nm(T: float, lambda_nm: np.ndarray) -> np.ndarray:
    """
    Planck spectral radiance B_λ(T) up to a constant factor across bands (relative shape only).
    λ input in nm. Output is relative (arbitrary units), sufficient for band weighting.
    """
    lam_m = np.asarray(lambda_nm, dtype=float) * 1e-9  # nm -> m
    lam_m = np.maximum(lam_m, 1e-20)
    # B_lambda = (2hc^2 / λ^5) * 1/(exp(hc/(λkT)) - 1)
    x = (_h * _c) / (lam_m * _kB * max(1e-12, float(T)))
    # Avoid overflow: for large x, exp(x) is huge; use np.exp with clipping
    x = np.clip(x, 1e-8, 1e3)
    denom = np.expm1(x)  # exp(x) - 1 with better small-x behavior
    B = (1.0 / (lam_m ** 5)) * (1.0 / (denom + 1e-30))
    # We omit the constant (2hc^2) because we only need relative shape
    B = np.clip(B, 0.0, np.inf)
    return B


def _normalize_spectrum_to_bands(B_lambda: np.ndarray, bands: SpectralBands) -> np.ndarray:
    """
    Normalize spectral radiance samples at band centers to produce band weights summing to 1
    using simple rectangle rule with band widths Δλ.
    """
    w = np.asarray(B_lambda, dtype=float) * np.asarray(bands.delta_lambda, dtype=float)
    wsum = float(np.sum(w)) + 1e-30
    return (w / wsum)


def blackbody_band_weights(T_eff: float, bands: SpectralBands) -> np.ndarray:
    """
    Compute normalized band weights (sum=1) for a blackbody at T_eff using band centers and widths.
    """
    Bc = _planck_lambda_nm(T_eff, bands.lambda_centers)
    return _normalize_spectrum_to_bands(Bc, bands)


def _rayleigh_band_factor(bands: SpectralBands) -> np.ndarray:
    mode = os.getenv("QD_ECO_TOA_TO_SURF_MODE", "simple").strip().lower()
    if mode != "rayleigh":
        return np.ones(bands.nbands, dtype=float)
    try:
        t0 = float(os.getenv("QD_ECO_RAYLEIGH_T0", "0.9"))
    except Exception:
        t0 = 0.9
    try:
        lref = float(os.getenv("QD_ECO_RAYLEIGH_LREF_NM", "550"))
    except Exception:
        lref = 550.0
    try:
        eta = float(os.getenv("QD_ECO_RAYLEIGH_ETA", "4.0"))
    except Exception:
        eta = 4.0
    return _rayleigh_weight(bands.lambda_centers, t0, lref, eta)


def dual_star_insolation_to_bands(insA: np.ndarray,
                                  insB: np.ndarray,
                                  bands: SpectralBands,
                                  *,
                                  T_eff_A: float | None = None,
                                  T_eff_B: float | None = None,
                                  j_A: float | None = None,
                                  j_B: float | None = None,
                                  L_ratio_A: float | None = None,
                                  M_ratio_A: float | None = None,
                                  L_ratio_B: float | None = None,
                                  M_ratio_B: float | None = None) -> np.ndarray:
    """
    Compute per-pixel band-averaged shortwave intensities [NB, nlat, nlon] for a dual-star system.

    Method:
      - Build per-star normalized spectra over bands from blackbody T_eff (or from L/M & j).
      - Combine at surface per pixel: S_b = (specA_b * insA + specB_b * insB) * Rayleigh(λ)
      - Normalize S_b over bands and scale by I_total = insA + insB. Nights (I_total≈0) stay zero.

    Args:
      insA, insB: 2D arrays [lat, lon] of per-star shortwave at surface/TOA proxy (W/m^2)
      bands: SpectralBands
      T_eff_*: optional known effective temperatures (K). If None, computed from (L_ratio, M_ratio, j).
      j_*: mass–radius exponent; defaults: env QD_STAR_A_J/QD_STAR_B_J (fallback 0.8)
      L_ratio_*, M_ratio_*: luminosity/mass ratio to Sun for each star; if None, read from constants.

    Returns:
      I_b: 3D array [NB, nlat, nlon] band intensities summing (≈) to (insA+insB).
    """
    # Lazy import to avoid cycles
    try:
        from pygcm import constants as const
    except Exception:
        const = None

    # Determine L/M ratios if not provided
    if L_ratio_A is None or M_ratio_A is None or L_ratio_B is None or M_ratio_B is None:
        if const is None:
            # conservative fallbacks ~ solar
            L_ratio_A = 1.0 if L_ratio_A is None else L_ratio_A
            M_ratio_A = 1.0 if M_ratio_A is None else M_ratio_A
            L_ratio_B = 1.0 if L_ratio_B is None else L_ratio_B
            M_ratio_B = 1.0 if M_ratio_B is None else M_ratio_B
        else:
            if L_ratio_A is None: L_ratio_A = float(getattr(const, "L_A", 1.0) / getattr(const, "L_SUN", 3.828e26))
            if M_ratio_A is None: M_ratio_A = float(getattr(const, "M_A", 1.0) / getattr(const, "M_SUN", 1.989e30))
            if L_ratio_B is None: L_ratio_B = float(getattr(const, "L_B", 1.0) / getattr(const, "L_SUN", 3.828e26))
            if M_ratio_B is None: M_ratio_B = float(getattr(const, "M_B", 1.0) / getattr(const, "M_SUN", 1.989e30))

    # j exponents
    if j_A is None:
        try:
            j_A = float(os.getenv("QD_STAR_A_J", "0.8"))
        except Exception:
            j_A = 0.8
    if j_B is None:
        try:
            j_B = float(os.getenv("QD_STAR_B_J", "0.8"))
        except Exception:
            j_B = 0.8

    # Effective temperatures
    if T_eff_A is None:
        # Allow explicit override
        env_TA = os.getenv("QD_STAR_A_TEFF_K")
        if env_TA:
            try:
                T_eff_A = float(env_TA)
            except Exception:
                T_eff_A = None
    if T_eff_B is None:
        env_TB = os.getenv("QD_STAR_B_TEFF_K")
        if env_TB:
            try:
                T_eff_B = float(env_TB)
            except Exception:
                T_eff_B = None
    if T_eff_A is None:
        T_eff_A = estimate_teff_from_LM(L_ratio_A, M_ratio_A, j=j_A, T_sun=_T_SUN)
    if T_eff_B is None:
        T_eff_B = estimate_teff_from_LM(L_ratio_B, M_ratio_B, j=j_B, T_sun=_T_SUN)

    # Per-star normalized spectra across bands
    specA = blackbody_band_weights(T_eff_A, bands)  # sum=1 over bands
    specB = blackbody_band_weights(T_eff_B, bands)

    # Atmospheric wavelength-dependent transmittance factor (Rayleigh/simple)
    T_ray = _rayleigh_band_factor(bands)  # shape [NB]
    T_ray = np.clip(T_ray, 0.0, np.inf)

    # Prepare shapes
    insA = np.asarray(insA, dtype=float)
    insB = np.asarray(insB, dtype=float)
    nlat, nlon = insA.shape
    NB = bands.nbands
    I_b = np.zeros((NB, nlat, nlon), dtype=float)

    # Combine per pixel:
    I_tot = insA + insB  # total shortwave (W/m^2)
    # Per-band raw signal at surface: (specA_b * insA + specB_b * insB) * T_ray_b
    # Then normalize over bands and scale by I_tot; nights remain zero.
    eps = 1e-12
    for b in range(NB):
        S_b = (specA[b] * insA + specB[b] * insB) * T_ray[b]
        I_b[b, :, :] = S_b

    # Normalize across bands and scale by I_tot
    S_sum = np.sum(I_b, axis=0)  # [nlat, nlon]
    # Avoid division by zero: where S_sum ~ 0 or I_tot ~ 0, set bands to zero.
    mask_pos = (S_sum > eps) & (I_tot > eps)
    if np.any(mask_pos):
        # Compute normalized fractions and scale by I_tot
        for b in range(NB):
            tmp = np.zeros_like(Sum := S_sum)
            tmp[mask_pos] = (I_b[b, :, :][mask_pos] / Sum[mask_pos]) * I_tot[mask_pos]
            I_b[b, :, :] = tmp
    else:
        I_b[:] = 0.0

    # Clean residual NaNs/Infs
    I_b = np.nan_to_num(I_b, nan=0.0, posinf=0.0, neginf=0.0)
    return I_b
