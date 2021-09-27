import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np

def path_to_mat(path, shape):
    """ Creates an rgb images based on a path of nodes.

    Arguments:
        path: a list of nodes.
        img: the image shape

    Returns:
        2D numpy image
    """
    mat = np.zeros(shape)
    mat = np.dstack([mat, mat, mat]) # rgb
    path = list(path)
    for i in range(len(path)):
        node = path[i]

        # if node.type == "Endpoint":
        #     mat[node.r, node.c, 0] = 255
        #     mat[node.r, node.c, 1] = 0
        #     mat[node.r, node.c, 2] = 0
        # else:
        mat[node.r, node.c,0] = 255
        mat[node.r, node.c,1] = 255
        mat[node.r, node.c,2] = 255

        if i == 0: 
            mat[node.r, node.c, 0] = 0
            mat[node.r, node.c, 1] = 255
            mat[node.r, node.c, 2] = 0
        elif i == len(path) - 1:
            mat[node.r, node.c, 0] = 0
            mat[node.r, node.c, 1] = 0
            mat[node.r, node.c, 2] = 255

    return mat

def binary_vector_to_bp_image(vector, width=50):
    """ Creates an image from an banding patterns.

    This basically expands the nubmer of rows.

    Arguments:
        vector: the banding pattern
        width: the number of rows the bp should be expanded to.

    Return:
        Banding pattern image
    """
    img = np.zeros((width, len(vector)))
    for i in range(len(vector)):
        val = vector[i]
        if val == 1:
            img[:, i] = 0
        elif val == 0:
            img[:, i] = 255
        elif val == -1:
            img[:, i] = 127
    
    return img.astype(int)

def binary_vector_comparison_img(vector1, vector2, width=50, buffer=5):
    """ Creates a matplotlib figure that comparing two list of banding patterns.

    Arguments:
        vector1: list of banding patterns.
        vecotor2: list of banding patterns.
        width: vertical length of each banding pattern
        buffer: space in between banding patterns.

    Returns:
        A matplotlib fig. Close the figure after saving etc.: fig.close()
    """
    n = int(np.sqrt(len(vector1)))

    imgs = []
    rows_img = width * 2 + buffer
    for i in range(len(vector1)):
        img1 = binary_vector_to_bp_image(vector1[i])
        img2 = binary_vector_to_bp_image(vector2[i])

        length = len(vector1[i])

        img_combined = np.zeros((rows_img, length))
        img_combined[0:width, :] = img1
        img_combined[width+buffer:, :] = img2
        imgs.append(img_combined)

    fig = plt.figure(figsize=(int(2*n), int(2*n)))
    grid = ImageGrid(fig, 111,
                     nrows_ncols=(n, n),
                     axes_pad=0.1)

    for ax, im in zip(grid, imgs):
        ax.imshow(im, cmap=plt.cm.gray)
        # ax[0].set_title('original')
        ax.axis('off')
    
    return fig
