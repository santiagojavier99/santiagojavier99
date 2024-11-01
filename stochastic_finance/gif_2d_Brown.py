import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import imageio

def create_plot(t, X, Y, x_lim, y_lim):
    plt.xlim(-x_lim, x_lim)
    plt.ylim(-y_lim, y_lim)
    plt.title(f'Frame {t}')  #f needed to read the string
    plt.scatter(X, Y)


    filename = f'frame_{t}.png'
    plt.savefig(filename)
    plt.close()

    return (filename)

def create_gif(filename_template, t_array, gif_filename='Brownian_2D.gif'):
    images=[]

    for t in t_array:
        current_filename=filename_template.format(t)
        images.append(imageio.imread(current_filename))

    imageio.mimsave(gif_filename, images, duratio=0.3)
    
