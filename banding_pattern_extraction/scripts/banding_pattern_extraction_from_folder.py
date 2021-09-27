import argparse
import os
import cv2 as cv
import csv
import matplotlib.pyplot as plt
from scipy.ndimage import binary_fill_holes
import traceback 

from .banding_pattern_extraction import get_banding_pattern

def folder_to_bp_csv(source_path, destination_path, extraction_size=None, identifier=None, csv_name=None, pixel_sampling=10, pixel_sigma=2, density_sigma=2, step_vector=1, chromsome_threshold=254, reject_multiple_blobs=False):
    """ Extracts the banding patterns from chromosomes in a folder to a csv file

    Arguments:
        source_path: path of the folder.
        destination_path: path to csv.
        exctation_size: size at which the banding pattern shall be extracted
        identifier: file identifier (only those will be considered). E.g. "23" only files with "23" in ther name will be extracted
        csv_name: name of the csv file.
        args**: see banding_pattern_extraction.py

    """
    
    file_list = os.listdir(source_path)
    banding_patterns = {}

    amount = len(file_list)

    # Extract all patterns, if possible, and save into dict
    for i in range(amount):
        file_name = file_list[i]

        if identifier != None and identifier not in file_name:
            continue

        if (i+1) % 100 == 0:
            print("Finished:", i / amount)
        file_path = os.path.join(source_path, file_name)
        img = cv.imread(file_path, 0)
        if extraction_size is not None:
            img = cv.resize(img, (extraction_size, extraction_size))

        try:
            bp = get_banding_pattern(img, pixel_sampling, pixel_sigma, density_sigma, step_vector, chromsome_threshold=chromsome_threshold, reject_multiple_blobs=reject_multiple_blobs)
            if not bp['error']:
                banding_patterns[file_name] = bp['binarized_banding_pattern']
            else:
                raise bp["error_message"]
        except Exception as e:
            print("Extraction of '{0}' failed, due to: {1}".format(file_name, e))
            traceback.print_exc()

    # Save dict with banding patterns
    if csv_name == None:
        csv_name = 'banding_patterns.csv'
    source_file = os.path.join(destination_path, csv_name)
    with open(source_file, 'w') as f:
        writer = csv.writer(f)

        # header
        writer.writerow(["file_name", "banding_pattern"])
        for file_name, bp in banding_patterns.items():
            bp_string = [str(x) for x in bp]
            bp_string = " ".join(bp_string)
            if identifier is not None:
                file_name = file_name.replace(identifier, '')
            writer.writerow([file_name, bp_string])