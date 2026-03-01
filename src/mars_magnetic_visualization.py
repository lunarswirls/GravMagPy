#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Imports:
import numpy as np
import pyshtools as pysh
import matplotlib.pyplot as plt

"""
@author lunarswirls
"""
# Average Mars radius
radius = pysh.constants.Mars.mean_radius  # meters

# Select ROI southwest of Medusa Fossae
lat_range = (-7, 3)
lon_range = (192, 202)

# SH coeffs from Langlais et al. 2019
coeffs_mag1 = pysh.datasets.Mars.Langlais2019()

# expand coeffs to maximum degree
B_grid1 = coeffs_mag1.expand()

fig1, ax1 = B_grid1.plot_total(cmap_limits=[0, 1250], show=False)
ax1.set_ylim([lat_range[0], lat_range[1]])
ax1.set_xlim([lon_range[0], lon_range[1]])
ax1.set_title('Medusa Fossae |B| - Langlais et al. 2019')
fig1.tight_layout()
fig1.show()

fig2, ax2 = B_grid1.plot_rad(cmap_limits=[-750, 750], show=False)
ax2.set_ylim([lat_range[0], lat_range[1]])
ax2.set_xlim([lon_range[0], lon_range[1]])
ax2.set_title('Medusa Fossae B_r - Langlais et al. 2019')
fig2.tight_layout()
fig2.show()

# SH coeffs from Morschhauser et al. 2014
coeffs_mag2 = pysh.datasets.Mars.Morschhauser2014()

# expand coeffs to maximum degree
B_grid2 = coeffs_mag2.expand()

fig3, ax3 = B_grid2.plot_total(cmap_limits=[0, 1250], show=False)
ax3.set_ylim([lat_range[0], lat_range[1]])
ax3.set_xlim([lon_range[0], lon_range[1]])
ax3.set_title('Medusa Fossae |B| - Morschhauser et al.')
fig3.tight_layout()
fig3.show()

fig4, ax4 = B_grid2.plot_rad(cmap_limits=[-750, 750], show=False)
ax4.set_ylim([lat_range[0], lat_range[1]])
ax4.set_xlim([lon_range[0], lon_range[1]])
ax4.set_title('Medusa Fossae B_r - Morschhauser et al.')
fig4.tight_layout()
fig4.show()

# SH coeffs from Wieczorek et al. 2024
coeffs_topo = pysh.datasets.Mars.MOLA_shape()

# expand coeffs to maximum degree of 5759
topo_grid = coeffs_topo.expand()
fig0, ax0 = topo_grid.plot(cmap_limits=[np.min(topo_grid.data)+0.0055*np.min(topo_grid.data), np.max(topo_grid.data)-0.0055*np.max(topo_grid.data)], show=False)
ax0.set_ylim([lat_range[0], lat_range[1]])
ax0.set_xlim([lon_range[0], lon_range[1]])
fig0.show()
