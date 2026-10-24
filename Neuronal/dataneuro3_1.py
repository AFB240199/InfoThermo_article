import numpy as np
import scipy.linalg as la
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from tqdm import tqdm
from matplotlib.colors import Normalize

# ============================================================
#  Dynamics + Jacobian
# ============================================================
def neuronal(t, x, A, F, B, R):
    phi = np.tanh(x)
    interaction = A @ phi
    return -B * x + F * phi + R * interaction


def jac_neuronal(t, x, A, F, B, R):
    N = len(x)
    I = np.eye(N)

    # sech^2(x) = 1 - tanh^2(x)
    sech2 = 1.0 - np.tanh(x) ** 2

    # Diagonal term
    J_diag = (-B + F * sech2) * I

    # Coupling term: d/dx_j [ R * sum_k A_ik tanh(x_k) ] = R * A_ij * sech^2(x_j)
    J_coupling = R * A * sech2[None, :]

    return J_diag + J_coupling


# ============================================================
#  Helpers
# ============================================================
def flip_signs_ensemble(A, z, p=0.26, base_seed=12345):
    """
    Create ensemble A_ens with same magnitudes/sparsity as A,
    but assign negative sign to a fraction p of the nonzero entries.
    Each realization uses a different seed (base_seed + i).
    """
    A = np.array(A, dtype=float, copy=True)
    nz = np.argwhere(A != 0.0)  # list of nonzero positions
    m = nz.shape[0]
    n_flip = int(np.round(p * m))

    A_ens = np.zeros((z, A.shape[0], A.shape[1]), dtype=float)
    seeds = base_seed + np.arange(z, dtype=int)

    for i, seed in enumerate(seeds):
        rng = np.random.default_rng(int(seed))
        idx = rng.choice(m, size=n_flip, replace=False)

        Ai = A.copy()
        rows = nz[idx, 0]
        cols = nz[idx, 1]

        # force chosen entries to be negative with same magnitude
        Ai[rows, cols] = -np.abs(Ai[rows, cols])

        A_ens[i] = Ai

    return A_ens, seeds


# --------------------------
# Laplacian (unsigned / absolute network)
# --------------------------
def build_laplacian_abs(A):
    """
    Option B (recommended baseline):
      A_abs = |A|
      D_ii  = sum_j A_abs_ij  (out-strength on absolute weights)
      L = D - A_abs
    """
    A_abs = np.abs(A)
    d = np.sum(A_abs, axis=1)
    D = np.diag(d)
    return D - A_abs


def entropy_from_U(U, eps=1e-12):
    """
    Von Neumann entropy of rho = U/Tr(U), with U symmetric real.
    """
    w = la.eigvalsh(U, check_finite=False)
    w = np.maximum(w, eps)
    w /= np.sum(w)
    return float(-np.sum(w * np.log(w)))


def evolve_entropy_neuronal(A_i, x_i, taus, dt, stride, F, B, R, eps=1e-12):
    """
    Given trajectory x_i[:,k], evolve U with time-dependent J(t) and return S(t_samp).
    """
    N = x_i.shape[0]
    taus_samp = taus[::stride]
    n_samp = len(taus_samp)
    S = np.full(n_samp, np.nan)

    U = np.eye(N)
    U /= np.trace(U)

    idx = 0
    for k in range(len(taus)):
        xk = x_i[:, k]
        t = taus[k]

        J = jac_neuronal(t, xk, A_i, F, B, R)

        # Liouville update
        Udot = J @ U + U @ J.T
        U += dt * Udot

        # symmetrize + real
        U = 0.5 * (U + U.T)
        U = U.real

        # normalize trace
        trU = float(np.trace(U).real)
        if abs(trU) < eps:
            trU = eps if trU >= 0 else -eps
        U /= trU

        if k % stride == 0:
            if not np.all(np.isfinite(U)):
                raise FloatingPointError("U contains NaN/Inf (neuronal)")
            S[idx] = entropy_from_U(U, eps=eps)
            idx += 1

    return S


def evolve_entropy_diffusive_from_L(L_abs, taus, dt, stride, R, eps=1e-12):
    """
    Diffusive baseline with Option B Laplacian:
      L_abs = D(|A|) - |A|
      J = -R * L_abs  (constant)
    """
    N = L_abs.shape[0]
    taus_samp = taus[::stride]
    n_samp = len(taus_samp)
    S = np.full(n_samp, np.nan)

    J = (-R * L_abs).real

    U = np.eye(N)
    U /= np.trace(U)

    idx = 0
    for k in range(len(taus)):
        Udot = J @ U + U @ J.T
        U += dt * Udot

        U = 0.5 * (U + U.T)
        U = U.real

        trU = float(np.trace(U).real)
        if abs(trU) < eps:
            trU = eps if trU >= 0 else -eps
        U /= trU

        if k % stride == 0:
            if not np.all(np.isfinite(U)):
                raise FloatingPointError("U contains NaN/Inf (diffusive)")
            S[idx] = entropy_from_U(U, eps=eps)
            idx += 1

    return S


def total_positive_time(S, dt_eff):
    """
    Total time where discrete derivative is positive:
    count of S[k] > S[k-1] times dt_eff.
    NaN-safe.
    """
    S = np.asarray(S)
    mask = np.isfinite(S[:-1]) & np.isfinite(S[1:])
    pos = (S[1:] > S[:-1]) & mask
    return float(np.sum(pos) * dt_eff)


# ============================================================
#  Main
# ============================================================
if __name__ == "__main__":

    # --------------------------
    # Load adjacency
    # --------------------------
    A = np.load("neuronal.npy")
    A = np.array(A, dtype=float)
    N = A.shape[0]

    # --------------------------
    # Parameters
    # --------------------------
    F = 0.0
    B = 1.0

    tf = 15.0
    n_steps = 1500
    dt = tf / n_steps
    taus = np.linspace(dt, tf, n_steps)

    eps = 1e-12
    stride = 15
    taus_samp = taus[::stride]
    n_samp = len(taus_samp)
    dt_eff = dt * stride

    # initial condition (fixed seed ONLY for x0)
    np.random.seed(1)
    x0 = np.random.normal(loc=0.02, scale=1.0 / N, size=N)

    # R values: at least 5 uniform in [0.05, 0.45]
    R_grid = np.linspace(0.05, 0.45, 5)
    n_R = len(R_grid)

    # Ensemble: sign-flip fraction p_neg on nonzero entries
    z = 30            # adjust as you like
    p_neg = 0.26
    A_ens, seeds_A = flip_signs_ensemble(A, z=z, p=p_neg, base_seed=12345)

    # ----------------------------------------------------------
    # Build ONE Laplacian on |A| (same for all A_i)
    # ----------------------------------------------------------
    L_abs = build_laplacian_abs(A)

    # Containers: (n_R, z, n_samp)
    S_neur = np.full((n_R, z, n_samp), np.nan)
    S_diff = np.full((n_R, z, n_samp), np.nan)

    # --------------------------
    # Compute DIFFUSIVE entropies ONCE per R, then broadcast over z
    # --------------------------
    for r_idx, R in enumerate(R_grid):
        try:
            Sd = evolve_entropy_diffusive_from_L(
                L_abs=L_abs, taus=taus, dt=dt, stride=stride, R=R, eps=eps
            )
            # replicate so downstream code can treat it like an ensemble
            S_diff[r_idx, :, :] = Sd[None, :]
        except Exception:
            # keep NaN for this R
            pass

    # --------------------------
    # Compute NEURONAL entropies per realization (needs trajectory)
    # --------------------------
    for r_idx, R in enumerate(R_grid):
        for i in tqdm(range(z), total=z, desc=f"R={R:.3f} | neuronal ensemble"):

            A_i = A_ens[i]

            try:
                # Trajectory (neuronal)
                sol = solve_ivp(
                    neuronal,
                    (dt, tf),
                    x0,
                    t_eval=taus,
                    args=(A_i, F, B, R),
                    # rtol=1e-6, atol=1e-9,
                )
                x_i = sol.y  # (N, T)

                # Entropy: neuronal (time-dependent J)
                S_neur[r_idx, i] = evolve_entropy_neuronal(
                    A_i=A_i, x_i=x_i, taus=taus, dt=dt, stride=stride,
                    F=F, B=B, R=R, eps=eps
                )

            except Exception:
                # keep NaN for this (r_idx, i)
                pass

    # --------------------------
    # Save (same structure as before, plus L_abs for reproducibility)
    # --------------------------
    np.savez_compressed(
        "entropy_neuronal_vs_diffusive_ensemble1.npz",
        taus=taus,
        taus_samp=taus_samp,
        stride=stride,
        dt=dt,
        dt_eff=dt_eff,
        eps=eps,
        z=z,
        p_neg=p_neg,
        seeds_A=seeds_A,
        F=F,
        B=B,
        R_grid=R_grid,
        S_neur=S_neur,
        S_diff=S_diff,
        # extra (harmless): makes clear what diffusive baseline was
        L_abs=L_abs,
    )

    # ============================================================
    #  Plot 1: total time with ΔS > 0 vs R (mean ± std, NaN-safe)
    # ============================================================
    posT_neur_mean = np.full(n_R, np.nan)
    posT_neur_std  = np.full(n_R, np.nan)
    posT_diff_mean = np.full(n_R, np.nan)
    posT_diff_std  = np.full(n_R, np.nan)

    for r_idx in range(n_R):
        Ts_n = []
        Ts_d = []

        for i in range(z):
            Sn = S_neur[r_idx, i]
            Sd = S_diff[r_idx, i]

            if np.all(np.isfinite(Sn)):
                Ts_n.append(total_positive_time(Sn, dt_eff))
            if np.all(np.isfinite(Sd)):
                Ts_d.append(total_positive_time(Sd, dt_eff))

        if len(Ts_n):
            posT_neur_mean[r_idx] = np.mean(Ts_n)
            posT_neur_std[r_idx]  = np.std(Ts_n, ddof=1) if len(Ts_n) > 1 else 0.0
        if len(Ts_d):
            posT_diff_mean[r_idx] = np.mean(Ts_d)
            posT_diff_std[r_idx]  = np.std(Ts_d, ddof=1) if len(Ts_d) > 1 else 0.0

    plt.figure(figsize=(10, 4))
    plt.errorbar(R_grid, posT_neur_mean, yerr=posT_neur_std, marker="o", lw=2, capsize=4, label="neuronal (J(t))")
    plt.errorbar(R_grid, posT_diff_mean, yerr=posT_diff_std, marker="o", lw=2, capsize=4, label="diffusive (J=-R L(|A|))")
    plt.xlabel("R")
    plt.ylabel("⟨ total time with ΔS>0 ⟩")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # ============================================================
    #  Plot 2: scatter per R
    #  Points = entropy values at different sampled times
    #  Color = time (early → late)
    #  Add diagonal y=x
    # ============================================================
    norm_t = Normalize(vmin=float(taus_samp[0]), vmax=float(taus_samp[-1]))

    for r_idx, R in enumerate(R_grid):

        Sn_mean_t = np.nanmean(S_neur[r_idx], axis=0)  # (n_samp,)
        Sd_mean_t = np.nanmean(S_diff[r_idx], axis=0)  # (n_samp,)

        mask = np.isfinite(Sn_mean_t) & np.isfinite(Sd_mean_t)
        x = Sn_mean_t[mask]
        y = Sd_mean_t[mask]
        c = taus_samp[mask]

        if len(x) == 0:
            continue

        Smin = float(np.min([x.min(), y.min()]))
        Smax = float(np.max([x.max(), y.max()]))
        pad = 0.05 * (Smax - Smin + 1e-12)

        plt.figure(figsize=(5.3, 5.0))
        sc = plt.scatter(x, y, c=c, s=28, alpha=0.9, norm=norm_t)

        xx = np.array([Smin - pad, Smax + pad])
        plt.plot(xx, xx, "--k", lw=1.1, alpha=0.7)

        plt.xlabel(r"$\langle S_{\mathrm{neur}}(t)\rangle_{\mathrm{ens}}$")
        plt.ylabel(r"$\langle S_{\mathrm{diff}}(t)\rangle_{\mathrm{ens}}$")
        plt.title(f"R = {R:.3f}")

        plt.axis("equal")
        plt.xlim(Smin - pad, Smax + pad)
        plt.ylim(Smin - pad, Smax + pad)
        plt.grid(True, alpha=0.3)

        cb = plt.colorbar(sc)
        cb.set_label("time t")

        plt.tight_layout()
        plt.show()
