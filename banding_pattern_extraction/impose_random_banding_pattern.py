import cv2 as cv
import numpy as np
import matplotlib.pyplot as plt
import argparse

from scripts.lib.utils import generate_random_banding_patterns
from scripts.chromosome_segmentation import get_segmented_chromosome
from scripts.impose_random_banding_pattern import impose_random_bp

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create banding pattern segmentation images of chromosome images.')
    parser.add_argument('-p', '--path', help='source path', required=True)
    args = parser.parse_args()

    img_path = args.path
    img = cv.imread(img_path, 0)
    fake_segmentation = impose_random_bp(img)
    real_segmentation = get_segmented_chromosome(img)
    
    fig, axes = plt.subplots(1, 3, figsize=(8, 8), sharex=True, sharey=True)
    ax = axes.ravel()

    ax[0].imshow(img, cmap=plt.cm.gray)
    ax[0].set_title("Banding Pattern")
    ax[0].axis('off')

    ax[1].imshow(real_segmentation, cmap=plt.cm.gray)
    ax[1].set_title("Real Band Segmentation")
    ax[1].axis('off')

    ax[2].imshow(fake_segmentation, cmap=plt.cm.gray)
    ax[2].set_title("Fake Band Segmentation")
    ax[2].axis('off')

    plt.show()