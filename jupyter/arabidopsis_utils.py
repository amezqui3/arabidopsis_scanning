import os
import math
import itertools

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy.ndimage as ndimage
import numpy.polynomial.polynomial as P

import tifffile as tf

import unionfind as UF
import MinimumBoundingBox as mbb

def persistence_with_UF(img, showfig=True, writefig=False, bname='intensity_histogram', dst='./', dpi=100):
    hist0,bins = np.histogram(img,bins=2**(img.dtype.itemsize*8),range=(0,2**(img.dtype.itemsize*8)))
    pers = sorted(UF.persistence(hist0[1:]),reverse=True)
    if showfig:
        plt.figure(figsize=(15,5))
        plt.plot(np.log(hist0[1:]+1), lw=3)
        plt.title(bname, fontsize=20);
        plt.xlabel('Intensity value')
        plt.ylabel('log(# voxels)')

        if writefig:
            filename = dst + bname + '.jpg'
            plt.savefig(filename, dpi=dpi, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight');
            plt.close()
    return pers

def clean_zeroes(img):
    d = img.ndim
    orig_size = img.size

    zeroes = []
    for ax in itertools.combinations(np.arange(d), d-1):
        row = np.any(img != 0, axis=ax)
        bar = np.nonzero(np.ediff1d(1-row))[0]
        row[bar] = True
        row[(bar+1)] = True
        zeroes.append(row)

    if d == 1:
        return img[zeroes[0]]

    mask = np.tensordot(zeroes[-1], zeroes[-2], axes=0).astype('bool')
    for i in range(3,d+1):
        mask = np.tensordot(mask, zeroes[-i], axes=0).astype('bool')

    shape = tuple([np.sum(zeroes[-i]) for i in range(1,len(zeroes)+1)])
    img = img[mask].reshape(shape)

    print(round(100-100*img.size/orig_size),'% reduction from input')

    return img, mask, shape

########################################################################
########################################################################

def get_individual_threshold(img, showfig=False):
    pers = persistence_with_UF(img, showfig)
    print(pers)

    if len(pers) >= 3:
        peaks = np.zeros(3, dtype=int)
        for i in range(len(peaks)):
            peaks[i] = pers[i][2]
        thr = int(np.average(peaks, weights=np.array([7,3,1])))
    if len(pers) == 2:
        thr = 0.1*(7*pers[0][2] + 3*pers[1][2])
    if len(pers) == 1:
        thr = pers[0][2]

    return thr

def plot4x4panel(img_list, ss, vmax=None, writefig=False, dst='./', bname='file'):

    rnum = 2
    cnum = len(img_list)//rnum
    fig, ax = plt.subplots(rnum,cnum,figsize=(5*cnum,5*rnum))
    print(rnum, cnum)

    if vmax is None:
        for idx in range(len(img_list)):
            i = idx//cnum
            j = idx%cnum
            ax[i,j].imshow(img_list[idx][ss], cmap='inferno', origin='lower');
    else:
        for idx in range(len(img_list)):
            i = idx//cnum
            j = idx%cnum
            ax[i,j].imshow(img_list[idx][ss], cmap='inferno', origin='lower', vmax = vmax[idx]);

    plt.tight_layout();

    if writefig:
        filename = pic_dst + bname + '_rind0.jpg'
        plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight');

def get_largest_element(comp):
    labels,num = ndimage.label(comp, structure=ndimage.generate_binary_structure(comp.ndim, 1))
    print(num,'components')
    hist,bins = np.histogram(labels, bins=num, range=(1,num+1))
    argsort_hist = np.argsort(hist)[::-1]
    print(np.sort(hist)[::-1][:20])

    j = 0
    i = argsort_hist[j]
    mask = labels==i+1
    box0 = comp.copy()
    box0[~mask] = 0
    box0[box0 > 0] = 1

    return box0

def fill_component(comp, x=True, y=True, z=True):
    rcomp = comp.copy()
    rcomp[rcomp > 0] = 1

    if x:
        for k in range(rcomp.shape[0]):
            rcomp[k,:,:] = ndimage.binary_fill_holes(rcomp[k,:,:])
        print('Closed X')
    if y:
        for k in range(rcomp.shape[1]):
            rcomp[:,k,:] = ndimage.binary_fill_holes(rcomp[:,k,:])
        print('Closed Y')
    if z:
        for k in range(rcomp.shape[2]):
            rcomp[:,:,k] = ndimage.binary_fill_holes(rcomp[:,:,k])
        print('Closed Z')

    return rcomp

def collapse_dimensions(img):
    snaps = []
    for i in range(img.ndim):
        snaps.append(np.sum(img, axis=i))
    return snaps

# same as above function, except take the max along each axis instead of the sum
def collapse_dimensions_max(img):
    snaps = []
    for i in range(img.ndim):
        snaps.append(np.max(img, axis=i))
    return snaps

# edit this function so that vmin and vmax can be specified while plotting
def plot_collapse_dimensions(snaps, bname='bname', tissue='tissue', display=False, writefig=False, dst='./', vmin = None, vmax = None):
    fig, ax = plt.subplots(1,len(snaps),figsize=(6*len(snaps),6))
    for i in range(len(snaps)):
        ax[i].imshow(snaps[i], cmap='inferno', origin='lower', vmin = vmin, vmax = vmax);
    plt.suptitle(bname + ' ' + tissue + ' collapse', fontsize=20);
    plt.tight_layout()

    if writefig:
        filename = dst + bname + '_' + '_'.join(tissue.split(' ')) + '.jpg'
        plt.savefig(filename, dpi=96, format='jpg', pil_kwargs={'optimize':True}, bbox_inches='tight');
        if not display:
            plt.close();
    return fig, ax
########################################################################
########################################################################


def tiff2coords(img, center=False):
    coords = np.nonzero(img)
    coords = np.vstack(coords).T
    if center:
        origin = -1*np.mean(coords, axis=0)
        coords = np.add(coords, origin)

    return coords

def plot_3Dprojections(seed, title='title', markersize=2, alpha=1, writefig=False, dst='./', dpi=150):
    axes = ['X','Y','Z']
    fig, ax = plt.subplots(1,3,figsize=(12,4))

    for i in range(3):
        proj = []
        for j in range(3):
            if j != i:
                proj.append(j)
        ax[i].plot(seed[:,proj[0]], seed[:,proj[1]], '.', ms=markersize, c='y', alpha=alpha)
        ax[i].set_xlabel(axes[proj[0]])
        ax[i].set_ylabel(axes[proj[1]])
        ax[i].set_title(axes[i] + ' Projection')
        ax[i].set_aspect('equal', 'datalim');

    fig.suptitle(title, y=0.95, fontsize=20)
    plt.tight_layout();

    if writefig:
        filename = '_'.join(title.split(' ')).lower()
        plt.savefig(dst + filename + '.png', dpi=dpi, format='png', bbox_inches='tight',
                    facecolor='white', transparent=False)
        plt.close();

########################################################################
########################################################################

def four_corners(boundary):
    coords = np.nonzero(boundary)
    coords = np.column_stack((coords[0], coords[1]))
    bbox = mbb.MinimumBoundingBox(coords)

    rect = np.array(list(bbox.corner_points))
    mincoord = np.min(rect, axis=0)
    maxcoord = np.max(rect, axis=0)

    tmost = rect[int(np.argwhere(rect[:,0] == maxcoord[0])[0]), :]
    rmost = rect[int(np.argwhere(rect[:,1] == maxcoord[1])[0]), :]

    bmost = rect[int(np.argwhere(rect[:,0] == mincoord[0])[0]), :]
    lmost = rect[int(np.argwhere(rect[:,1] == mincoord[1])[0]), :]

    return np.array([tmost, rmost, bmost, lmost])

def four_borderlines(corners):
    coefs = np.empty((len(corners), corners.ndim))

    for i in range(len(corners)):
        j = (i+1)%len(corners)

        line = np.vstack((corners[i,:], corners[j,:]))
        poly_fit = P.Polynomial.fit(line[:,1] , line[:,0], deg=1, full=False )
        coefs[i,:] = poly_fit.convert().coef

    return coefs

def four_midpoints(corners):
    midpoints = np.empty((len(corners), corners.ndim))

    for i in range(len(corners)):
        j = (i+1)%len(corners)

        midpoints[i,:] = (corners[i,:] + corners[j,:])/2
    return midpoints

def midline_splits(midpoints):
    splits = np.empty((int(len(midpoints)/2), midpoints.ndim))
    for i in range(len(splits)):
        j = (i+2)%len(midpoints)

        line = np.vstack((midpoints[i,:], midpoints[j,:]))
        poly_fit = P.Polynomial.fit(line[:,1] , line[:,0], deg=1, full=False )
        splits[i,:] = poly_fit.convert().coef

    return splits

def four_centers(corners, midpoints):
    center = np.mean(corners,axis=0)
    centers = np.empty((len(corners), corners.ndim))

    for i in range(len(corners)):
        j = (i-1)%len(corners)
        pts = np.vstack((center, corners[i,:], midpoints[i,:], midpoints[j,:]))
        centers[i,:] = np.mean(pts, axis=0)

    return centers

def normalize_density(img, adjust_by):
    resol = 2**(img.dtype.itemsize*8)
    npz = np.arange(resol, dtype=img.dtype)

    for i in range(len(npz)):
        aux = round(adjust_by[0] + adjust_by[1]*npz[i])
        if aux < resol and aux > 0:
            npz[i] = int(aux)
        elif aux >= resol:
            npz[i] = resol - 1
        else:
            npz[i] = 0

    with np.nditer(img, flags=['external_loop'], op_flags=['readwrite']) as it:
        for x in it:
            x[...] = npz[x]

    return img
