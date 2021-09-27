import numpy as np
from scipy.ndimage import gaussian_filter1d

def interpolate_ends(r, c, blobs):
    """ Linerally interpolate the ends of a skeleton.

    r and c are going to be extended as along as one index outside of the given segmentation.

    Arguments:
        r: row indices.
        c: column indices.
        blobs: chromosome binary segmentation

    Returns:
        Interpolate tuple of r and c.
    """
    r_interpolated, c_interpolated = _interpolate(r, c, blobs, end=False)
    r_interpolated, c_interpolated = _interpolate(r_interpolated, c_interpolated, blobs, end=True)
    
    return r_interpolated, c_interpolated

def _interpolate(r_vec, c_vec, blobs, end=True):
    """ Linerally interpolates given list of indices on one side

    Arguments:
        r_vec: row indices.
        c_vec: column indices.
        blobs: chromosome binary segmentation
        end: optional, bool, whether the end should be interpolated

    Returns:
        The onsided interpolation r_vec and c_vec
    """
    if end:
        r1, r2 = r_vec[-2:]
        c1, c2 = c_vec[-2:]
    else:
        r2, r1 = r_vec[:2]
        c2, c1 = c_vec[:2]

    step_c = c2 - c1
    step_r = r2 - r1

    # TODO: check for identical value

    # norm = np.sqrt(step_c**2 + step_r**2)
    # step_c /= norm
    # step_r /= norm

    new_r = []
    new_c = []

    count = 1
    skeleton_not_reached = True
    while skeleton_not_reached:
        
        current_r = int(np.round(r2 + step_r * count))
        current_c = int(np.round(c2 + step_c * count))
        try:
            skeleton_not_reached = blobs[current_r, current_c]
        except Exception: # out of bounds
            skeleton_not_reached = False

        new_r.append(current_r)
        new_c.append(current_c)

        if count > 5000:
            raise Exception("Too many iterations")
            # TODO: Catch this error early on
            # Comment: should not occurr anymore, but will keep this, just in case

        count += 1

    new_r = np.array(new_r)
    new_c = np.array(new_c)

    if end:
        r_vec = np.append(r_vec, new_r)
        c_vec = np.append(c_vec, new_c)
    
    else:
        r_vec = np.append(np.flip(new_r), r_vec)
        c_vec = np.append(np.flip(new_c), c_vec)
    
    return r_vec, c_vec

def smoothen(r, c, sigma):
    """ Smooths a given list of row and column indices

    Arguments:
        r_vec: row indices.
        c_vec: column indices.
        sigma: sigma of Gaussian kernel

    Returns:
        A tuple of smoothed r and c
    """
    r = gaussian_filter1d(r, sigma)
    c = gaussian_filter1d(c, sigma)

    # Remove same values (can happen after smoothing)
    r_new = [r[0]]
    c_new = [c[0]]
    for i in range(1, len(r)):
        if r_new[-1] != r[i] or c_new[-1] != c[i]:
            r_new.append(r[i])
            c_new.append(c[i])

    r = np.array(r_new)
    c = np.array(c_new)

    return r, c

def subsample(r, c, sampling):
    """ Subsamples row and column indices.

    Arguments:
        r: row indices.
        c: column indices.
        sampling: Keep every x-th pixel (if possible), otherwise keep ends
    
    Returns:
        Subsampled r and c
    """
    amount_values = len(r)
    sample_value = int(amount_values / sampling)

    if sample_value > 1: # We have at least two values
        indx_mask = np.linspace(0, len(r) - 1, sample_value).astype(int)
        r = r[indx_mask]
        c = c[indx_mask]
    else: # simply return first and last value
        r = np.array([r[0], r[-1]])
        c = np.array([c[0], c[-1]])

    return r, c