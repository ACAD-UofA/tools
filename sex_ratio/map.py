#!/usr/bin/env python

from __future__ import print_function
import sys
import collections

import numpy as np
import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def parse_csv(fn):
    with open(fn) as f:
        header = next(f).rstrip().split("\t")
        for i, line in enumerate(f,2):
            line = line.rstrip()
            fields = line.split("\t")
            mfields = {}
            for h, f in zip(header, fields):
                mfields[h] = f

            try:
                sex = mfields["sex"]

                try:
                    lat = float(mfields["lat"])
                    lon = float(mfields["lon"])
                except ValueError:
                    continue

                try:
                    age = float(mfields["age"])
                except ValueError:
                    age = None
            except:
                print("Exception at {}:{}".format(fn, i), file=sys.stderr)
                raise
            yield sex, lat, lon, age

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: {} omap.pdf bison.csv brownbears.csv".format(sys.argv[0]),
                file=sys.stderr)
        exit(1)

    omap_fn = sys.argv[1]
    ifiles = sys.argv[2:]

    scale = 0.025
    h, w = 210.0, 297.0
    #w, h = 4, 3
    hcentroid = 160
    jitter = 5
    proj = ccrs.PlateCarree(central_longitude=hcentroid)
    geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))

    fig1 = plt.figure(figsize=(scale*w,scale*h))
    gs1 = gridspec.GridSpec(3, 2, width_ratios=[5,1], height_ratios=[5,5,1])
    ax1 = fig1.add_subplot(gs1[0,:], projection=proj)
    ax2 = fig1.add_subplot(gs1[1,:], projection=proj)
    ax3 = fig1.add_subplot(gs1[2,0])
    ax4 = fig1.add_subplot(gs1[2,1])

    for ifile_fn, title, ax in zip(ifiles, ("Bison", "Brown bears"), (ax1, ax2)):

        lat = collections.defaultdict(list)
        lon = collections.defaultdict(list)
        age = collections.defaultdict(list)
        for _sex, _lat, _lon, _age in parse_csv(ifile_fn):
            theta = np.random.uniform(0, 2*np.pi)
            r = np.random.uniform(0, jitter/2.0)
            xjitter, yjitter = r*np.cos(theta), r*np.sin(theta)
            _lon += 360-hcentroid + xjitter
            _lat += yjitter
            lat[_sex].append(_lat)
            lon[_sex].append(_lon)
            age[_sex].append(_age)

        #ax.set_global()
        ax.set_extent((0,330, 0, 90), geodetic)

        ax.add_feature(cfeature.LAND)
        #ax.add_feature(cfeature.OCEAN)
        ax.add_feature(cfeature.COASTLINE, lw=0.1, zorder=-1)
        #ax.add_feature(cfeature.BORDERS, linestyle=':', lw=0.5)

        #ax.gridlines()

        for sex, m in zip("MF", "oX"):
            cb = ax.scatter(lon[sex], lat[sex],
                    marker=m,
                    #facecolor=fc,
                    edgecolors="black",
                    c=age[sex],
                    #alpha=0.8,
                    s=15,
                    lw=0.5,
                    cmap="plasma",
                    norm=matplotlib.colors.LogNorm(vmin=100, vmax=50000),
                    label=sex,
                    )

        ax.set_title(title)
        #ax.legend(framealpha=1.0)

    leg = [Line2D([0], [0], marker='o', markeredgecolor='black', label='Male',
                          markerfacecolor='white', markersize=10, lw=0),
           Line2D([0], [0], marker='X', markeredgecolor='black', label='Female',
                          markerfacecolor='white', markersize=10, lw=0)]

    clb = plt.colorbar(cb, cax=ax3, orientation='horizontal')
    clb.ax.set_title("Radiocarbon age (years)")

    ax4.legend(handles=leg, loc="center", frameon=False)
    ax4.axis("off")

    fig1.tight_layout()
    fig1.savefig(omap_fn)
