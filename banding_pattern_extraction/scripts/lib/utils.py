import numpy as np
import noise

def clip_bp(bp, bp_max_length):
    """ Clips a banding pattern symmetrically.

    Arguments:
        bp: The banding pattern.
        bp_max_length: new length

    Returns:
        The clipped banding pattern
    """
    half_diff = (len(bp)- bp_max_length) / 2

    if bp_max_length < 0:
        raise ValueError("No clipping needed, bp small enough")

    left_pad = int(np.floor(half_diff))
    right_pad = int(np.ceil(half_diff))

    end = int(len(bp) - right_pad)
    
    return bp[left_pad:end]

def pad_bp(bp, bp_max_length, pad_value):
    """ Pads a banding pattern symmetrically.

    Arguments:
        bp: The banding pattern.
        bp_max_length: new length.
        pad_values: numerical pad value

    Returns:
        The padded banding pattern
    """
    padded = np.ones(bp_max_length) * pad_value
    half_diff = (bp_max_length - len(bp)) / 2

    left_pad = int(np.floor(half_diff))
    right_pad = int(np.ceil(half_diff))

    end = int(len(padded) - right_pad)
    padded[left_pad:end] = bp

    return padded.astype(int)

def one_hot_encode(bp, num_labels=3):
    """ On hot encodes a banding pattern

    Arguments:
        bp: The banding pattern.
        num_labels: number of classes.


    Returns:
        A one hot encoded representation of the banding patterns
    """
    one_hot = np.zeros((bp.size, num_labels))
    one_hot[np.arange(bp.size), bp] = 1

    return one_hot

def generate_random_banding_patterns(batch_size, max_bp_length, class_statistics, scale=10, octaves=1, lacunarity=3, persistence=0.7):
    """ Generates a random banding pattern based in perlin noise

    Arguments:
        batch_size: amoung ot banding patterns to create
        max_bp_length: the maximum length of a banding pattern
        class statistics: A list of tuples that each contain the mean and variance of a class

    Returns:
        - The padded perlin noise banding patterns as a numpy array
        - The one hot encoded perlin noise banding patterns as a numpy array
        - List of length of each banding pattern
    """

    bps_input = np.zeros((batch_size, max_bp_length))
    bps_one_hot = np.zeros((batch_size, max_bp_length, 3))
    lengths = []
    
    amount_modes = len(class_statistics)
    for j in range(batch_size):

        mode = np.random.randint(0, amount_modes)
        mean_length = class_statistics[mode][0]
        std_length = class_statistics[mode][1]
        # Randomly sample the length of the chromosome based on passed mean and std
        rnd_length = int(np.random.normal(mean_length, std_length))
        rnd_length = np.max([np.min([rnd_length, max_bp_length]), 0])
        lengths.append(rnd_length)

        # Random starting point for perlin noise
        base = np.random.randint(np.iinfo(np.int8).min, np.iinfo(np.int8).max)

        # Fill vector with perlin values
        bp = np.zeros(rnd_length)
        for i in range(rnd_length):
            bp[i] = noise.pnoise1((i + 1) / scale, octaves=octaves, lacunarity=lacunarity, persistence=persistence, base=base, repeat=np.iinfo(np.int16).max)

        # Perlin noise between -1 and 1
        bp = bp > 0 
        # Pad with -1 at the sides, save as input for network
        bp_input = pad_bp(bp, max_bp_length, -1)
        bps_input[j,:] = bp_input

        # Create categorical one hot equivalent, where the -1 corresponds to class 2, for loss calculation
        bp_categorical = bp_input.copy()
        bp_categorical[np.where(bp_categorical==-1)] = 2
        one_hot_bp = one_hot_encode(bp_categorical)
        bps_one_hot[j,:,:] = one_hot_bp

    return bps_input, bps_one_hot, lengths