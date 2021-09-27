from skimage.morphology import medial_axis, skeletonize, skeletonize_3d
from multiprocessing import Process, Queue, current_process, Manager, Pool
from time import time
import queue
import sys
import matplotlib.pyplot as plt
import cv2 as cv
import numpy as np
from scipy.ndimage import binary_fill_holes
import argparse
import traceback

from .lib.Graph import Graph
from .lib.path_preprocessing import interpolate_ends, smoothen, subsample
from .lib.banding_pattern_utils import *

def get_banding_pattern(img, pixel_sampling=8, pixel_sigma=2, density_sigma=2, step_vector=1, chromsome_threshold=254, size=None, reject_multiple_blobs=False, pickle_conform_results=False):
    """ Extracs the banding pattern of a stained chromosome image.

    Arguments:
        img: the chromosome image.
        pixel_sampling: optional, Sub sampling rate, i.e. keep every x-th pixel from the interpolated skeleton.
        pixel_sigma: CURRENTLY NOT IN USE, optinal, Sigma for Gaussian filter that is applied on the skeleton pixels pixel prior sub sampling.
        density_sigma: optional, Sigma for Gaussian filter of the density profiles.
        chromsome_threshold: optional, the grayscale segmentation threshold. 
        size: optional, size of the extracted banding pattern. Will be rescaled! Does not change the size of the image.
        reject_multiple_blobs: optional, stops the algorithm early if multiply blobs have been detected.
        pickle_conform_results: optional, only return serializable results.


    Returns:
        A dictionary with many intermediate results, look into the command line interface for more info.
        Final banding pattern can be accesses by: dicts["binarized_banding_pattern"]
    """

    try:

        if img is None:
            raise ValueError("Empty image")

        # Create skeleton
        if chromsome_threshold == None:
            chromsome_threshold = np.median(img) 
        
        blobs = img < chromsome_threshold
        blobs = binary_fill_holes(blobs)
        skeleton = skeletonize_3d(blobs) # lee method

        # Create graph of the skeleton, and get its nodes
        graph = Graph(skeleton)

        # If there is only a single pixel in the skeleton, add another one on top of it.
        # This will make the algorithm create a vertical medial axis. We asssume
        # that the chromsomes are aligned along the vertical axiss.
        if graph.max_cluster_length < 2:
            node = list(graph.nodes.values())[0]
            skeleton[node.r + 1 , node.c] = 255
            graph = Graph(skeleton)

        # Get number of blobs detected (clusters in the graph)
        num_blobs = graph.amount_clusters

        # Get shortest paths from each endpoint to each other endpoint in the biggest cluster
        paths = graph.endpoint_paths

        # Some error handling
        if num_blobs == 0:
            raise ValueError("No path detected. Either empty image or circular structure.")

        if num_blobs > 1 and reject_multiple_blobs:
            raise ValueError("Multiple blobs detected, early rejection.")
        
        if graph.endpoints == 0:
            raise ValueError("Circular structure detected. Cluster only has one endpoint.")

        # Merge single paths, to remove branches
        longest_path = graph.get_longest_merged_path(paths)

        # Extract the row and columns from the nodes
        r, c = graph.nodes_to_numpy(longest_path)

        # TODO: Maybe include smoothing again
        # ## Smoothen values
        # r, c = smoothen(r, c, pixel_sigma)

        ## Subsample the vertices, use linspace, so we dont drop the last value
        r, c = subsample(r, c, pixel_sampling)

        # if len(r) > 3:
        #     r = r[1:len(r)-1]
        #     c = c[1:len(c)-1]

        # Interpolate ends
        r_interpolated, c_interpolated = interpolate_ends(r, c, blobs)

        # Sample pixels across
        r_sampled, c_sampled, banding_points, banding_pattern = sample(r_interpolated, c_interpolated, blobs, img, res=step_vector)

        # Flip banding pattern if input was upside down
        if r[0] > r[-1]:
            banding_pattern = np.flip(banding_pattern)

        # Filter banding pattern
        banding_pattern_filtered, banding_pattern_smooth = banding_pattern_filter(banding_pattern, density_sigma)

        # Binarize Banding Pattern
        binarized_banding_pattern = binarize_banding_pattern(banding_pattern_filtered)

        # Potentially resize bp
        if size is not None:
            binarized_banding_pattern = resize_banding_pattern(binarized_banding_pattern, size)

        if pickle_conform_results:
            results = {
                'binarized_banding_pattern': binarized_banding_pattern,
                'num_blobs': num_blobs,
                'error': False,
                'error_message': ''
            }

        else:
            results = {
                'binarized_banding_pattern': binarized_banding_pattern,
                'banding_pattern_filtered': banding_pattern_filtered,
                'banding_pattern_smooth': banding_pattern_smooth,
                'banding_points': banding_points,
                'banding_pattern': banding_pattern,
                'r': r,
                'c': c,
                'r_interpolated': r_interpolated,
                'c_interpolated': c_interpolated,
                'r_sampled': r_sampled,
                'c_sampled': c_sampled,
                'paths': paths,
                'longest_path': longest_path,
                'skeleton': skeleton,
                'blobs': blobs,
                'num_blobs': num_blobs,
                'error': False,
                'error_message': '',
                'stack_trace': ''
            }

        return results

    except Exception as e:

        if hasattr(e, 'message'):
            message = e.message
        else:
            message = "No error message"

        results = {
            'error': True,
            'error_message': message,
            'stack_trace':  traceback.format_exc()
        }

        return results


def get_banding_pattern_multi_process(imgs, workers, pixel_sampling=5, pixel_sigma=2, density_sigma=2, step_vector=1, chromsome_threshold=None, size=None, reject_multiple_blobs=False):
    """ Extracts multiple banding patterns with multiple processes

    Arguments:
        img: numpy 3D array, where the first dimension corresponds to the sample
        workers: number of processes
        args**: get_banding_pattern()

    Returns:
        List of results, where each entry contains intermediate values and the final binarized banding pattern.
        See get_banding_pattern()
    """
    
    amount_imgs = imgs.shape[0]
    manager = Manager()
    pending_imgs = Queue()
    results = manager.Queue()
    processes = []
    params = {
        "pixel_sampling":pixel_sampling, 
        "pixel_sigma":pixel_sigma, 
        "density_sigma":density_sigma, 
        "step_vector":step_vector, 
        "chromsome_threshold":chromsome_threshold, 
        "size":size,
        "reject_multiple_blobs":reject_multiple_blobs,
        "pickle_conform_results": True, # Very important, some of the variables are not serialiasable
    }

    for i in range(amount_imgs):
        pending_imgs.put((i, imgs[i], params))

    # creating processes
    for w in range(workers):
        p = Process(target=__queue_get_banding_pattern, args=(pending_imgs, results))
        processes.append(p)
        p.start()

    # completing process
    for p in processes:
        p.join()

    # print the output
    output_results = [0] * amount_imgs
    while not results.empty():
        num, result = results.get()
        output_results[num] = result

    return output_results

def __queue_get_banding_pattern(pending_imgs, results):
    """ Helper function for handling a queue of processes extracting banding patterns

    Arguments:
        pending_imgs: list of images.
        results: list where results are saved.
        
    """
    while True:
        try:
            num, img, params = pending_imgs.get_nowait()
        except queue.Empty:
            break
        else:
            result = get_banding_pattern(img, **params)
            results.put((num, result))
    return True