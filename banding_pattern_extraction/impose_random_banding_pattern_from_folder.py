import argparse
import os
import cv2 as cv
import traceback 

from scripts.impose_random_banding_pattern import impose_random_bp

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extracts the banding patterns from chromosome images from a folder.')
    parser.add_argument('-s', '--source_path', help='source path', required=True)
    parser.add_argument('-d', '--destination_path', help='destination path', required=True)

    args = parser.parse_args()
    source_path = args.source_path
    destination_path = args.destination_path

    file_list = os.listdir(source_path)
    banding_patterns = {}

    amount = len(file_list)

    # Extract all patterns, if possible, and save into dict
    for i in range(amount):
        file_name = file_list[i]
        if (i+1) % 100 == 0:
            print("Finished:", i / amount)
        file_path = os.path.join(source_path, file_name)
        save_path = os.path.join(destination_path, file_name)
        
        img = cv.imread(file_path, 0)

        try:
            bp_segmentation = impose_random_bp(img)
            cv.imwrite(save_path, bp_segmentation)
            

        except Exception as e:
            print("Extraction of '{0}' failed, due to: {1}".format(file_name, e))
            traceback.print_exc()

