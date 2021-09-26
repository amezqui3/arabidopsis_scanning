import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import tifffile as tf
import os
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_path',help='input tiff file or directory')
    parser.add_argument('min_size',help='minimum w*h*d',nargs='?',default=0,type=int)
    args = parser.parse_args()

    files = []
    if os.path.isdir(args.in_path):
        for r,_,fnames in os.walk(args.in_path):
            for fname in fnames:
                files.append((r+'/'+fname))
    elif os.path.isfile(args.in_path):
        files.append(args.in_path)

    resol,stepsize,C = 32,16,0

    print('Found {} file{}'.format(len(files),'' if len(files)==1 else 's'))
    for f in files:

        img = tf.imread(f)
        d,h,w = img.shape
        print(f,w,h,d,args.min_size,img.nbytes//1e6)
        #img = tf.transpose_axes(img, 'YXZ', 'XYZ')
        d,h,w = img.shape
        aux = 1
        cmap = matplotlib.cm.get_cmap('inferno',stepsize)
        bounds = list(range(0,resol, resol//stepsize))
        bounds.append(resol)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        if C:
            npz = np.arange(resol).astype('uint8')
            z11 = np.array([3.17108944e-04, 9.92360110e-01, 1.61116136e+00])
            for i in range(len(npz)):
                aux = round(z11[0]*npz[i]*npz[i]+z11[1]*npz[i]+z11[2])
                if aux < 256 and aux > 0:
                    npz[i] = int(aux)
                elif aux > 255:
                    npz[i] = 255
                else:
                    npz[i] = 0

            with np.nditer(img, flags=['external_loop'], op_flags=['readwrite']) as it:
                for x in it:
                    x[...] = npz[x]

        print(f,w,h,d,args.min_size)

        if d*h*w>=args.min_size:
            tf.imshow(img,cmap = cmap, norm = norm, title=f,origin='lower')
            # tf.imshow(img,title=f,cmap=plt.cm.magma)
            plt.show()

        C = 1

if __name__ == '__main__':
    main()
