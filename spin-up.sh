#!/usr/bin/env bash
# Qingdai GCM rapid spin-up (two-phase SOP)
# Phase 1: shallow mixed layer (fast equilibration), save restart
# Phase 2: restore physical MLD, warm-restart to equilibrium
set -euo pipefail

# -------- Configurable knobs (override via env) --------
PHASE1_YEARS="${PHASE1_YEARS:-150}"         # years in Phase 1 (shallow MLD)
PHASE2_YEARS="${PHASE2_YEARS:-120}"         # years in Phase 2 (physical MLD)
MLD_P1="${MLD_P1:-5}"                       # m, shallow mixed-layer depth for Phase 1
MLD_P2="${MLD_P2:-50}"                      # m, physical mixed-layer depth for Phase 2
RESTART_P1="${RESTART_P1:-restart_phase1.nc}"
RESTART_EQ="${RESTART_EQ:-restart_equilibrium.nc}"

# Optional: pick latest topography NetCDF if available (unless QD_TOPO_NC is pre-set)
if [[ -z "${QD_TOPO_NC:-}" ]]; then
  latest_nc="$(ls -t data/*.nc 2>/dev/null | head -n1 || true)"
  if [[ -n "$latest_nc" ]]; then
    export QD_TOPO_NC="$latest_nc"
    export QD_USE_TOPO_ALBEDO="${QD_USE_TOPO_ALBEDO:-1}"
  fi
fi

# -------- Common recommended settings (can be overridden by env) --------
export QD_ENERGY_W="${QD_ENERGY_W:-1}"
export QD_ENERGY_DIAG="${QD_ENERGY_DIAG:-1}"
export QD_HUMIDITY_DIAG="${QD_HUMIDITY_DIAG:-1}"
export QD_WATER_DIAG="${QD_WATER_DIAG:-1}"

export QD_USE_OCEAN="${QD_USE_OCEAN:-1}"
export QD_OCEAN_USE_QNET="${QD_OCEAN_USE_QNET:-1}"
export QD_USE_SEAICE="${QD_USE_SEAICE:-1}"

# Dynamics filtering (P010)
export QD_FILTER_TYPE="${QD_FILTER_TYPE:-combo}"
export QD_SIGMA4="${QD_SIGMA4:-0.02}"
export QD_SHAPIRO_EVERY="${QD_SHAPIRO_EVERY:-6}"
export QD_SHAPIRO_N="${QD_SHAPIRO_N:-2}"

# Plot cadence (reduce I/O for long runs)
export QD_PLOT_EVERY_DAYS="${QD_PLOT_EVERY_DAYS:-20}"

# -------- Phase 1: Rapid Equilibration (shallow MLD) --------
echo "== Phase 1: Rapid Equilibration =="
export QD_INIT_BANDED="${QD_INIT_BANDED:-1}"          # apply banded initial Ts
export QD_INIT_T_EQ="${QD_INIT_T_EQ:-295.0}"
export QD_INIT_T_POLE="${QD_INIT_T_POLE:-265.0}"

export QD_MLD_M="${MLD_P1}"
export QD_TOTAL_YEARS="${PHASE1_YEARS}"
unset QD_RESTART_IN
export QD_RESTART_OUT="${RESTART_P1}"

echo "Phase1 config: YEARS=${QD_TOTAL_YEARS}, MLD=${QD_MLD_M} m, RESTART_OUT=${QD_RESTART_OUT}"
python3 -m scripts.run_simulation

# -------- Phase 2: Full Physics Adjustment (physical MLD) --------
echo "== Phase 2: Full Physics Adjustment =="
export QD_MLD_M="${MLD_P2}"
export QD_TOTAL_YEARS="${PHASE2_YEARS}"
export QD_RESTART_IN="${RESTART_P1}"
export QD_RESTART_OUT="${RESTART_EQ}"

echo "Phase2 config: YEARS=${QD_TOTAL_YEARS}, MLD=${QD_MLD_M} m, RESTART_IN=${QD_RESTART_IN}, RESTART_OUT=${QD_RESTART_OUT}"
python3 -m scripts.run_simulation

echo "Spin-up complete."
echo "Phase-1 restart: ${RESTART_P1}"
echo "Equilibrium restart: ${RESTART_EQ}"
