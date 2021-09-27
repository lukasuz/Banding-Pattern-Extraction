import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.interpolate import interp1d

from .bresenham import bresenham_pixel_summation

# Rotation vector for 90 degrees
R = np.array([
    [np.cos(np.pi/2), -np.sin(np.pi/2)],
    [np.sin(np.pi/2), np.cos(np.pi/2)]
])

def in_bounds(img, r, c):
    """ Checks whether a row or column value is in bounds in the image

    Arguments:
        img: the image.
        r: row index.
        c: column index.

    Returns:
        bool
    """
    return r < img.shape[0] and r >= 0 and c < img.shape[1] and c >=0

def sample(r_points, c_points, blobs, img, res=1, max_length=50):
    """ Samples the grayscale values for each perpendicular line all given point in a chromosome.

    Arguments:
        r_points: array of row indices
        c_points: array of column indices
        blobs: the segmentation of the chromosome image
        res: optional, "resolution" in what frequencies we step along the medial axis. Anything else than 1 
            breaks the code for now.
        max_length: optional, max_length * 2 is the maximum length of the perpendicular line that is going
            to be sampled. If you have very big images, increase this.
        
    Returns:
        Tuple of four objects: 
            1. list of row indices on the medial axis that were used for sampling
            2. list of column indices on the medial axis that were used for sampling
            3. List of list, where each sublist contains the image indices that were used for sampling
                one entry in the denity profile.
            4. The raw density profile (mean grayscale level across each perpendicular line)

        new_points[:,0], new_points[:, 1], banding_points, banding_pattern
    """
    points = np.array([r_points, c_points]).T
    vectors = points[1:,:] - points[:-1,:]
    norm = np.sqrt(vectors[:,0]**2 + vectors[:,1]**2)

    norm_vectors = vectors / np.expand_dims(norm, axis=-1)
    step_vectors = norm_vectors * res
    perpendicular_vectors = norm_vectors @ R
    middle_points = norm / 2
    
    new_points = []
    banding_pattern = []
    banding_points = [[0,0]]
    prev_buffer = 0
    for i in range(points.shape[0] - 1): # iterate over each point
        fixed_norm = norm[i] - prev_buffer
        samples_per_vector = fixed_norm // res
        middle_point = middle_points[i]

        steps = np.arange(prev_buffer, samples_per_vector + res, res)
        for step in steps:
            new_point = points[i] + step_vectors[i] * step

            new_point_pixel = np.round(new_point).astype(int)
            try:
                if not blobs[new_point_pixel[0], new_point_pixel[1]] or not in_bounds(img, new_point_pixel[0], new_point_pixel[1]):
                    continue
            except Exception: # Out of bounds
                break
            

            new_points.append(new_point)

            current_vector = perpendicular_vectors[i]
            if step < middle_point:

                try:
                    neighbour_vector = perpendicular_vectors[i-1]
                    middle_point_prev = middle_points[i-1]
                except Exception:
                    neighbour_vector = current_vector
                    middle_point_prev = middle_point
                
                dist = middle_point + middle_point_prev 
                c1 = (middle_point - step) / dist
                c2 = 1 - c1

            else:
                
                try:
                    neighbour_vector = perpendicular_vectors[i+1]
                    middle_point_next = middle_points[i+1]
                except Exception:
                    neighbour_vector = current_vector
                    middle_point_next = middle_point

                dist = middle_point + middle_point_next
                c1 = (step - middle_point) / dist
                c2 =  1 - c1

            banding_vector = current_vector * c2 + neighbour_vector * c1
            # banding_vector = current_vector

            end_point_1 = new_point + banding_vector * max_length
            end_point_2 = new_point - banding_vector * max_length

            banding_indx_1, banding_val_1 = bresenham_pixel_summation(new_point, end_point_1, blobs, img)
            banding_indx_2, banding_val_2 = bresenham_pixel_summation(new_point, end_point_2, blobs, img)
            banding_pattern.append(np.mean([banding_val_1[1:] + banding_val_2]))

            banding_points.append(list(reversed(banding_indx_1[1:])) + banding_indx_2)
   

        remainder_per_vector = np.remainder(fixed_norm, res)

        if i + 1 == vectors.shape[0]:
            continue

        if remainder_per_vector != 0.0:
            if np.all(norm_vectors[i] == norm_vectors[i+1]):
                prev_buffer = res - remainder_per_vector
            else:
                C = np.arccos(np.clip(np.dot(norm_vectors[i], norm_vectors[i+1]), -1.0, 1.0))
                C = np.pi - C
                sin_C = np.sin(C)
                x = remainder_per_vector * sin_C
                A = np.arcsin(x / res)
                B = np.pi - C - A
                y = res * np.sin(B)
                prev_buffer = y / sin_C
        else:
            prev_buffer = 0

    new_points = np.array(new_points)
    # banding_points = np.array(banding_points)
    # new_vec = new_points[1:,:] - new_points[:-1,:]
    # new_norm = np.sqrt(new_vec[:,0]**2 + new_vec[:,1]**2)
    # print(new_norm)
    return new_points[:,0], new_points[:, 1], banding_points, banding_pattern

def banding_pattern_filter(banding_pattern, sigma=2):
    """ Applies a gaussian filter and the non linear filter to the raw density profile
    
    Arguments:
        banding_pattern: the raw banding pattern, i.e. the density profile
        sigma: optional, the sigma of the gaussian kernel

    Returns:
        The filtered density profile.

    """
    
    # f_b = median_filter(banding_pattern, sigma)
    f_b = gaussian_filter1d(banding_pattern, sigma)
    

    prev_i_b = np.copy(f_b)
    i_b = []
    R = 2
    same_pattern = False
    counter = 0
    while not same_pattern:
        counter += 1
        for i in range(0, len(f_b)):
            
            crrnt = prev_i_b[i] 
            if i == 0:
                prv = crrnt
            else:
                prv = prev_i_b[i-1]
            
            if i == len(f_b) - 1:
                nxt = crrnt
            else:
                nxt = prev_i_b[i+1]

            neighbourhood = [prv, crrnt, nxt]
            dif_min = crrnt - np.min(neighbourhood)
            dif_max = np.max(neighbourhood) - crrnt

            if dif_max <= dif_min:
                i_b.append(crrnt + dif_max / R)
            else:
                i_b.append(crrnt - dif_min / R)

        R = max(R - 1, 1)
        same_pattern = np.all(prev_i_b == i_b)
        prev_i_b = i_b
        i_b = []

    return prev_i_b, f_b


def binarize_banding_pattern(banding_pattern, black_tag=1, white_tag=0):
    """ Binarizes a density profile
    
    Arguments:
        banding_pattern: the filtered density profile.
        black_tag: the value that valleys (black bands) should be taged with.
        white_tag: the value that peaks (white bands) should be taged with.

    Returns:
        The binarized banding pattern.
    
    """
    b1 = mark_diffs(banding_pattern, black_tag=black_tag, white_tag=white_tag)
    b2 = np.flip(mark_diffs(np.flip(banding_pattern), black_tag=black_tag, white_tag=white_tag))
    b3 = np.abs(np.array(b1)-np.array(b2))

    clusters = cluster_1D(b3)
    final_bp = _retag_saddle_point(clusters, b2)

    return final_bp

def resize_banding_pattern(banding_pattern, new_length):
    """ Resizes a banding pattern by the means of linear interpolation
    
    Argument:
        banding pattern: Any type of banding pattern.
        new_length: the wished length.

    Returns:
        The resized banding pattern.

    """
    length = len(banding_pattern)

    x = np.linspace(0, length, length)
    f = interp1d(x, banding_pattern)
    x_new = np.linspace(0, length, new_length)
    resized_banding_pattern = np.round(f(x_new))

    return resized_banding_pattern

def cluster_1D(arr):
    """ Clusters a 1D binary vector

    Arguments:
        arr: a binary vector

    Returns:
        A dictinary containing the index starting point of a cluster and its length.
    """
    clusters = {}
    i = 0
    while i < len(arr):
        if arr[i] == 1:
            j = i
            while j < len(arr) and arr[j] == 1:
                j += 1
            clusters[i] = j - i
            i = j
        else:
            i += 1
    return clusters

def mark_diffs(banding_pattern, white_tag = 1, black_tag = 0):
    """ Helper function. Marks the saddle points in a filtered density profile.
    
    """
    binarized_bp = []
    first_band = True
    current_tag = 0
    for i in range(len(banding_pattern)):

        if i == 0:
            diff = 0
        else:
            diff = banding_pattern[i-1] - banding_pattern[i]
        
        if diff > 0: # positive diff, we have a vally
            current_tag = black_tag
            # Tag all prev. values as white band, if we are in a valley 
            if first_band:
                binarized_bp = [white_tag] * len(binarized_bp)
                first_band = False

        elif diff < 0: # negative diff, we have a peak
            current_tag = white_tag
            # tag all prev. values as black band if we are in a peal
            if first_band:
                binarized_bp = [black_tag] * len(binarized_bp)
                first_band = False
        else: 
            pass # grad == 0, dont change tag
        binarized_bp.append(current_tag)
    
    return binarized_bp


def _retag_saddle_point(clusters, source):
    """ Biniarization helper. Retags saddle points.
    """
    b4 = np.copy(source)
    for indx, amount in clusters.items():
        half = amount // 2
        b4[indx:indx+half] = source[indx-1]
        b4[indx+half:indx+amount] = source[indx+amount]
    return b4