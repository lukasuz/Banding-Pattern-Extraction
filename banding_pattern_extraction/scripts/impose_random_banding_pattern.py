import cv2 as cv
import numpy as np
import matplotlib.pyplot as plt
import argparse

from .lib.utils import generate_random_banding_patterns
from .banding_pattern_extraction import get_banding_pattern
from .chromosome_segmentation import nearest_nonzero_idx, get_segmented_chromosome

def impose_random_bp(img, extraction_size=None, fake_bp=None, **args):
    """ Imposes a random Perlin noise banding pattern onto a chromosome image shape.

    Arguments:
        img: the chromosome image.
        extraction_size: optional, size, at which the chromsome BP shall be extracted and imposed (128x128)
            works well, assumes square image.
        fake_bp: optional, banding pattern that should be imposed, otherwise a random one
            is generated.
        args**: see banding_pattern_extraction.py
    
    Returns:
        The segmentation of the imposed banding pattern onto the chromosome.
    """



    original_size = img.shape
    if extraction_size is not None:
        img = cv.resize(img, (extraction_size, extraction_size))

    chromosome_segmentation = np.zeros_like(img).astype(float)
    counter = np.zeros_like(img)

    results = get_banding_pattern(img, **args)
    real_bp = results['binarized_banding_pattern']
    bp_points = results['banding_points']
    blob = results['blobs']

    length_bp = len(real_bp)

    if fake_bp is None:
        fake_bp = np.squeeze(generate_random_banding_patterns(1, length_bp, [[length_bp, 0]])[0])
    
    for i in range(length_bp):
        banding_pattern_line = np.array(bp_points[i])
        try:
            r = banding_pattern_line[:,1]
            c = banding_pattern_line[:,0]
            chromosome_segmentation[c, r] += (1 - fake_bp[i])
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
        final_segmentation = cv.resize(final_segmentation, (original_size[1], original_size[0]))

    return final_segmentation
