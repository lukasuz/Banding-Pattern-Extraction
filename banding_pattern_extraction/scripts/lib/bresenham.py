import numpy as np

def bresenham_pixel_summation(p1, p2, blobs, img):
    """ Adaption of bresenham algorithm that retrieves the grayscale values across a line

    Adaption of https://github.com/encukou/bresenham/blob/master/bresenham.py

    Arguments:
        p1: point one
        p2: point two
        blobs: chromosome binary segmentation image
        img: the chromosome image
    
    Returns:
        A tuple:
            1. The indices
            2. List of grayscale values
    """

    x0, y0 = np.round(p1).astype(int)
    x1, y1 = np.round(p2).astype(int)
    dx = x1 - x0
    dy = y1 - y0

    xsign = 1 if dx > 0 else -1
    ysign = 1 if dy > 0 else -1

    dx = abs(dx)
    dy = abs(dy)

    if dx > dy:
        xx, xy, yx, yy = xsign, 0, 0, ysign
    else:
        dx, dy = dy, dx
        xx, xy, yx, yy = 0, ysign, xsign, 0

    D = 2*dy - dx
    y = 0

    grayscale_vals = []
    indx = []
    for x in range(dx + 1):
        x_temp = x0 + x*xx + y*yx
        y_temp = y0 + x*xy + y*yy

        try:
            if not blobs[x_temp, y_temp]:
                break
            else:
                grayscale_vals.append(img[x_temp, y_temp])
                indx.append([x_temp, y_temp])
        except Exception: # out of bounds
            break
        if D >= 0:
            y += 1
            D -= 2*dx
        D += 2*dy

    return indx, grayscale_vals