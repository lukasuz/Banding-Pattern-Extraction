import cv2 as cv
import numpy as np
import matplotlib.pyplot as plt
import argparse

from .banding_pattern_extraction import get_banding_pattern

def get_segmented_chromosome(img, extraction_size=None, **args):
    """ Creates a chromosome segmentation mask based on the extracted chromosome

    Arguments:
        img: the chromosome image.
        extraction_size: optional, size, at which the chromsome BP shall be extracted and imposed (128x128)
            works well, assumes square image.
        **args: see banding_pattern_extraction.py
    Returns:
        The banding pattern chromosome segmentation mask.
    """

    original_size = img.shape
    if extraction_size is not None:
        img = cv.resize(img, (extraction_size, extraction_size))

    chromosome_segmentation = np.zeros_like(img)
    counter = np.zeros_like(img)

    results = get_banding_pattern(img, **args)
    bp = results['binarized_banding_pattern']
    bp_points = results['banding_points']
    blob = results['blobs']

    for i in range(len(bp)):
        banding_pattern_line = np.array(bp_points[i])
        try:
            r = banding_pattern_line[:,1]
            c = banding_pattern_line[:,0]
            chromosome_segmentation[c, r] += (1 - bp[i])
            counter[c, r] += 1

        except: # Empty line, simply ignore
            pass
    
    divison_mask = np.where(counter > 0)

    chromosome_segmentation[divison_mask] = chromosome_segmentation[divison_mask] / counter[divison_mask]
    chromosome_segmentation = np.round(chromosome_segmentation)
    chromosome_segmentation = chromosome_segmentation * 127.5 # white bands should be gray

    final_segmentation = np.ones_like(img) * 255 # background should be white
    final_segmentation[divison_mask] = chromosome_segmentation[divison_mask]

    # Now where have to fill in values that were not sampled
    xs, ys = np.where(blob == 1)

    for i in range(len(xs)):
        x = xs[i]
        y = ys[i]

        if final_segmentation[x, y] < 255:
            continue

        x_neighbour, y_neighbour = nearest_nonzero_idx(blob, x, y)
        final_segmentation[x, y] = final_segmentation[x_neighbour, y_neighbour]
    
    if extraction_size is not None:
        final_segmentation = cv.resize(final_segmentation, original_size)

    return final_segmentation


def nearest_nonzero_idx(a,x,y):
    """ https://stackoverflow.com/questions/43306291/find-the-nearest-nonzero-element-and-corresponding-index-in-a-2d-numpy-array """
    idx = np.argwhere(a)

    # If (x,y) itself is also non-zero, we want to avoid those, so delete that
    # But, if we are sure that (x,y) won't be non-zero, skip the next step
    idx = idx[~(idx == [x,y]).all(1)]

    return idx[((idx - [x,y])**2).sum(1).argmin()]

