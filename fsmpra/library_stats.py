"""Library-level statistics shared by review and site-methylation stages."""
import numpy as np
from scipy import stats


def split_half(a, b):
    """Split-half correlation, Spearman-Brown reliability and signal/noise SDs (pp)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    r = float(np.corrcoef(a, b)[0, 1])
    cov = float(np.cov(a, b)[0, 1])
    return dict(n=len(a), pearson_halves=r, spearman_halves=float(stats.spearmanr(a, b)[0]),
                reliability_full=2 * r / (1 + r), signal_sd_pp=float(np.sqrt(max(cov, 0.))),
                noise_sd_full_pp=float(np.sqrt(np.var(a - b, ddof=1) / 4)))


def position_means(pos, values):
    pos, values = np.asarray(pos), np.asarray(values, float)
    keep = np.isfinite(values)
    upos = np.unique(pos[keep])
    return upos, np.array([values[keep & (pos == p)].mean() for p in upos])


def circular_shift(pos, chromatin, activity, margin=10):
    """Exact circular position-shift test of position-level correlation.

    Position means keep the three substitutions at a site together; rotating one
    profile against the other preserves each profile's spatial autocorrelation.
    Every rotation with offset in [margin, L-margin] is enumerated, so the
    smallest attainable p is 1/(number of rotations + 1)."""
    pos, chromatin, activity = np.asarray(pos), np.asarray(chromatin, float), np.asarray(activity, float)
    keep = np.isfinite(chromatin) & np.isfinite(activity)
    pc, xc = position_means(pos[keep], chromatin[keep])
    pa, ya = position_means(pos[keep], activity[keep])
    if not np.array_equal(pc, pa):
        raise ValueError('Position sets differ')
    observed = float(np.corrcoef(xc, ya)[0, 1])
    shifts = range(margin, len(xc) - margin + 1)
    null = np.array([np.corrcoef(np.roll(xc, s), ya)[0, 1] for s in shifts])
    p = (1 + np.sum(np.abs(null) >= abs(observed))) / (len(null) + 1)
    return dict(positions=len(xc), position_pearson=observed, rotations=len(null),
                null_abs_max=float(np.abs(null).max()), p_circular_shift=float(p))
