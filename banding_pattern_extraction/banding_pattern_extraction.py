from time import time
import matplotlib.pyplot as plt
import cv2 as cv
import numpy as np
import argparse
import traceback

from scripts.banding_pattern_extraction import get_banding_pattern
from scripts.lib.banding_pattern_utils import *

if __name__ == "__main__":
    from scripts.lib.visualisation_utils import *

    parser = argparse.ArgumentParser(description='Extract the banding pattern from a chromosome image.')
    parser.add_argument('-p', '--path', help='input image', required=True)
    parser.add_argument('--pixel_sampling', help='Sub sampling rate, i.e. keep every x-th pixel from the interpolated skeleton', type=int, nargs='?', default=8)
    parser.add_argument('--pixel_sigma', help='Sigma for Gaussian filter that is applied on the skeleton pixels pixel prior sub sampling', type=int, nargs='?', default=2)
    parser.add_argument('--density_sigma', help='Sigma for Gaussian filter of the density profiles', type=int, nargs='?', default=2)
    parser.add_argument('--threshold', help='Segmentation threshold', type=int, nargs='?', default=254)
    parser.add_argument('--size', help='Size of the  banding pattern, extraction will be resized accordingly', type=int, nargs='?', default=None)

    args = parser.parse_args()
    img_path = args.path
    pixel_sampling = args.pixel_sampling
    pixel_sigma = args.pixel_sigma
    density_sigma = args.density_sigma
    chromsome_threshold = args.threshold
    size = args.size
    step_vector = 1

    t1 = time()

    ## Load img
    img = cv.imread(img_path, 0)

    # A lot of results for visualisation purposes
    try:
        t1 = time()
        results = get_banding_pattern(img, pixel_sampling, pixel_sigma, density_sigma, step_vector, chromsome_threshold=254, size=size)
        print("Time: ", time() - t1)
        
        if results["error"]:
            print(results["error_message"])
            print(results["stack_trace"])

    except Exception:
        traceback.print_exc()
        exit()
        
    
    binarized_banding_pattern = results['binarized_banding_pattern']
    banding_pattern_filtered = results['banding_pattern_filtered']
    banding_pattern_smooth = results['banding_pattern_smooth']
    banding_points = results['banding_points']
    banding_pattern = results['banding_pattern']
    r = results['r']
    c = results['c']
    r_interpolated = results['r_interpolated']
    c_interpolated = results['c_interpolated'] 
    r_sampled = results['r_sampled']
    c_sampled = results['c_sampled']
    paths = results['paths']
    longest_path = results['longest_path']
    skeleton = results['skeleton']
    blobs = results['blobs']
    num_blobs = results['num_blobs']

    print("Length of banding pattern:", len(binarized_banding_pattern))
    print("Number of clusters: ", num_blobs)

    t2 = time()
    print("Seconds taken:", t2 - t1)

    # Plot unfiltered density profile
    fig = plt.figure()
    x = np.arange(len(banding_pattern))
    plt.bar(x, banding_pattern)
    plt.title("Density Profile")

    # Plot filtered density profile
    fig = plt.figure()
    x = np.arange(len(banding_pattern))
    plt.bar(x, banding_pattern_filtered)
    plt.title("Filtered Density Profile")

    # Plot Banding Pattern
    fig = plt.figure()
    bp_img = binary_vector_to_bp_image(binarized_banding_pattern)
    plt.imshow(bp_img, cmap=plt.cm.gray)
    plt.title("Banding Pattern")

    # Combined 
    fig, axes = plt.subplots(1, 4, figsize=(8, 8), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(img, cmap=plt.cm.gray)
    ax[0].set_title('original')
    ax[0].axis('off')

    ax[1].imshow(skeleton, cmap=plt.cm.gray)
    ax[1].set_title("Skeleton")
    ax[1].axis('off')

    path_img = path_to_mat(longest_path, skeleton.shape)
    ax[2].imshow(path_img)
    ax[2].contour(blobs, cmap=plt.cm.gray)
    ax[2].scatter(c, r, s=5, c='r')
    ax[2].set_title('Longest Path')
    ax[2].axis('off')

    ax[3].imshow(blobs, cmap=plt.cm.gray)
    ax[3].plot(c_interpolated, r_interpolated)
    ax[3].scatter(c_sampled, r_sampled, s=3, c='r')
    ax[3].scatter(c_interpolated, r_interpolated, s=3, c='g')
    for i in range(len(banding_points)):
        banding_pattern_line = np.array(banding_points[i])
        try:
            ax[3].plot(banding_pattern_line[:,1], banding_pattern_line[:,0])
        except Exception:
            pass
    ax[3].scatter(c, r, s=20, c='b')
    ax[3].set_title('Sampled')
    ax[3].axis('off')

    # Plot paths etc.
    amount = int(np.ceil(len(paths) / 2))
    fig, axes = plt.subplots(amount, 5, figsize=(8, 8), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(img, cmap=plt.cm.gray)
    ax[0].set_title('original')
    ax[0].axis('off')

    ax[1].imshow(blobs, cmap=plt.cm.gray)
    ax[1].plot(c_interpolated, r_interpolated)
    ax[1].scatter(c_sampled, r_sampled, s=3, c='r')
    ax[1].scatter(c_interpolated, r_interpolated, s=3, c='g')
    for i in range(len(banding_points)):
        banding_pattern_line = np.array(banding_points[i])
        try:
            ax[1].plot(banding_pattern_line[:,1], banding_pattern_line[:,0])
        except Exception:
            pass
    ax[1].scatter(c, r, s=20, c='b')
    ax[1].set_title('Interpolated')
    ax[1].axis('off')

    path_img = path_to_mat(longest_path, skeleton.shape)
    ax[2].imshow(path_img)
    ax[2].contour(blobs, cmap=plt.cm.gray)
    ax[2].scatter(c, r, s=5, c='r')
    ax[2].set_title('Longest Path')
    ax[2].axis('off')

    for i in range(len(paths)):
        path = paths[i]
        path_img = path_to_mat(path, skeleton.shape)

        ax[i+3].imshow(path_img, cmap='magma')
        ax[i+3].contour(blobs, [0.5], colors='w')
        ax[i+3].set_title('path: ' + str(i))
        ax[i+3].axis('off')

    fig.tight_layout()
    plt.show()