import argparse
import os
import csv
import traceback 

from scripts.banding_pattern_extraction import get_banding_pattern

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extracts the banding patterns from chromosome images from a folder.')
    parser.add_argument('-s', '--source_path', help='source path', required=True)
    parser.add_argument('-d', '--destination_path', help='destination path', required=True)
    parser.add_argument('--pixel_sampling', help='Sub sampling rate, i.e. keep every x-th pixel from the interpolated skeleton', type=int, nargs='?', default=5)
    parser.add_argument('--pixel_sigma', help='Sigma for Gaussian filter that is applied on the skeleton pixels pixel prior sub sampling', type=int, nargs='?', default=2)
    parser.add_argument('--density_sigma', help='Sigma for Gaussian filter of the density profiles', type=int, nargs='?', default=2)
    parser.add_argument('--threshold', help='Segmentation threshold', type=int, nargs='?', default=254)

    args = parser.parse_args()
    source_path = args.source_path
    destination_path = args.destination_path
    pixel_sampling = args.pixel_sampling
    pixel_sigma = args.pixel_sigma
    density_sigma = args.density_sigma
    chromsome_threshold = args.threshold
    step_vector = 1

    file_list = os.listdir(source_path)
    banding_patterns = {}

    amount = len(file_list)

    # Extract all patterns, if possible, and save into dict
    for i in range(amount):
        file_name = file_list[i]
        if (i+1) % 100 == 0:
            print("Finished:", i / amount)
        file_path = os.path.join(source_path, file_name)
        img = cv.imread(file_path, 0)

        try:
            bp = get_banding_pattern(img, pixel_sampling, pixel_sigma, density_sigma, step_vector, chromsome_threshold=chromsome_threshold, reject_multiple_blobs=False)
            if not bp['error']:
                banding_patterns[file_name] = bp['binarized_banding_pattern']
            else:
                raise bp["error_message"]
        except Exception as e:
            print("Extraction of '{0}' failed, due to: {1}".format(file_name, e))
            traceback.print_exc()

    # Save dict with banding patterns
    source_file = os.path.join(destination_path, 'banding_patterns.csv')
    with open(source_file, 'w') as f:
        writer = csv.writer(f)

        # header
        writer.writerow(["file_name", "banding_pattern"])
        for file_name, bp in banding_patterns.items():
            bp_string = [str(x) for x in bp]
            bp_string = " ".join(bp_string)
            writer.writerow([file_name, bp_string])