import argparse
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import cv2 as cv
import sys

from .scripts.banding_pattern_extraction import get_banding_pattern

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extracts Statistics of the banding patterns from chromosome images in a give folder.')
    parser.add_argument('-s', '--source_path', help='source path', required=True)
    parser.add_argument('-d', '--destination_path', help='destination path', required=True)

    args = parser.parse_args()
    source_path = args.source_path

    data_path = args.source_path
    saving_path = args.destination_path

    os.makedirs(saving_path, exist_ok=True)

    # Statistics By Type
    if(os.path.isdir(data_path)):

        lens = []
        i = 0
        f = open(saving_path + "/types.csv",'x')
        for type in range(1,24):
            lens = []
            max_bp  = 0
            min_bp  = 0
            mean_bp = 0
            total   = 0
            for filename in os.listdir(data_path):
        
                if filename.endswith(".png") and filename.startswith(str(type)+'_'): #by type
                #if filename.endswith(".png"):

                    print('Processing: ' + filename)
                    img_path = data_path + '/' + filename
                    img = cv.imread(img_path, 0)
                    results = get_banding_pattern(img, chromsome_threshold=254)
                    #print(results['binarized_banding_pattern'])
                    
                    if(not results['error']):
                        lens.append(len(results['binarized_banding_pattern']))
                        i = i+1

            max_bp  = max(lens)
            min_bp  = min(lens)
            mean_bp = np.mean(lens)
            total   = len(lens)
            print('TYPE '+ str(type))
            print('-----------------')
            print('max = ' + str(max_bp)) #Max is 256 for type 12
            print('min = ' + str(min_bp))
            print('mean = ' + str(mean_bp))
            print('total = ' + str(total))
            # save to file
            f.write('TYPE '+ str(type) + '\n')
            f.write('----------\n')
            f.write('max = ' + str(max_bp) + '\n')
            f.write('min = ' + str(min_bp) + '\n')
            f.write('mean = ' + str(mean_bp) + '\n')
            f.write('total = ' + str(total) + '\n')
            f.write('Lengths vector:' + '\n')
            f.write(str(lens) + '\n')
            f.write('\n')
            
        f.close()


    # Statistics by All images
    if(os.path.isdir(data_path)):
        
        lens = []
        i = 0
        f = open(saving_path + "/all.csv",'x')

        max_bp  = 0
        min_bp  = 0
        mean_bp = 0
        total   = 0

        for filename in os.listdir(data_path):
        
            #if filename.endswith(".png") and filename.startswith(str(type)+'_'): #by type
            if filename.endswith(".png"):

                print('Processing: ' + filename)
                img_path = data_path + '/' + filename
                img = cv.imread(img_path, 0)
                results = get_banding_pattern(img, chromsome_threshold=254)
                #print(results['binarized_banding_pattern'])
                    
                if(not results['error']):
                    lens.append(len(results['binarized_banding_pattern']))
                    i = i+1

        max_bp  = max(lens)
        min_bp  = min(lens)
        mean_bp = np.mean(lens)
        total   = len(lens)
        print('All ')
        print('-----------------')
        print('max = ' + str(max_bp)) #Max is 256 for type 12
        print('min = ' + str(min_bp))
        print('mean = ' + str(mean_bp))
        print('total = ' + str(total))
        # save to file
        f.write('max = ' + str(max_bp) + '\n')
        f.write('min = ' + str(min_bp) + '\n')
        f.write('mean = ' + str(mean_bp) + '\n')
        f.write('total = ' + str(total) + '\n')
        f.write('Lengths vector:' + '\n')
        f.write(str(lens) + '\n')
        f.write('\n')

        f.close()

        print('Done!')

    else:
        print('Wrong path to data!')