import cv2 as cv
import matplotlib.pyplot as plt
import argparse

from scripts.chromosome_segmentation import get_segmented_chromosome

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create banding pattern segmentation images of chromosome images.')
    parser.add_argument('-p', '--path', help='source path', required=True)
    parser.add_argument('--extraction_size', default=None, type=int, help='scaling factor for extraction of banding factor')
    args = parser.parse_args()

    img_path = args.path
    img = cv.imread(img_path, 0)
    segmentation = get_segmented_chromosome(img, extraction_size=args.extraction_size)
    
    fig = plt.figure()
    plt.imshow(img, cmap=plt.cm.gray)
    plt.title("Banding Pattern")

    fig = plt.figure()
    plt.imshow(segmentation, cmap=plt.cm.gray)
    plt.title("Banding Segmentation")
    plt.show()