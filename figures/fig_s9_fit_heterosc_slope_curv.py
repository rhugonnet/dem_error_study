"""Plotting of Figure S9: 1D fit for heteroscedasticity with slope or curvature for the Mont-Blanc case study"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Open estimates
fn_bin = '/home/atom/ongoing/work_stderr_dem/case_study_montblanc/df_heteroscedas_slope_curv.csv'

df_bin = pd.read_csv(fn_bin)

fig = plt.figure(figsize=(6,7))

# First, fit with slope
grid = plt.GridSpec(22, 23, wspace=0.1, hspace=0.1)

ax = fig.add_subplot(grid[:10, :])

df_slp = df_bin[np.logical_and(df_bin.nd==1, np.isfinite(df_bin.slope_mid))]
df_slp = df_slp
def slope_f(x,a,b):
    return min(df_slp.nmad)+a*np.exp(-b*x/np.pi)

cof, _ = curve_fit(slope_f, df_slp.slope_mid.values[:-1], df_slp.nmad.values[:-1], method='trf')

x = np.linspace(0, 90)

ax.scatter(df_slp.slope_mid, df_slp.nmad, marker='x')
ax.plot(x, slope_f(x, *cof), linestyle='dashed', color='black')
ax.set_xlabel('Slope (degrees)')
ax.set_ylabel('Dispersion (1$\sigma$) of\nelevation differences (m)')

ax.text(0.5, 0.9, '$f(x) = a + b \cdot e^{-cx}$' + '\nwith a = {:.1f}, b = {:.1f} and c = {:.1f}'.format(min(df_slp.nmad),cof[0], cof[1]*np.pi), transform=ax.transAxes, ha='center', va='top')
ax.text(0.025, 0.95, 'a', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

# Second, fit with curvature
ax = fig.add_subplot(grid[13:, :])

df_maxc = df_bin[np.logical_and(df_bin.nd==1, np.isfinite(df_bin.maxc_mid))]

p = np.polyfit(df_maxc.maxc_mid.values, df_maxc.nmad.values, 1)
x = np.linspace(0, 1.2*max(df_maxc.maxc_mid))

ax.scatter(df_maxc.maxc_mid, df_maxc.nmad, marker='x')
ax.plot(x, np.polyval(p,x), linestyle='dashed', color='black')
ax.set_xlabel('Maximum curvature categories (10$^{2}$ m$^{-1}$)')
ax.set_ylabel('Dispersion (1$\sigma$) of\nelevation differences (m)')
ax.text(0.5, 0.9, '$f(x) = a \cdot x + b$' +'\nwith a = {:.1f}, b = {:.1f}'.format(p[0], p[1]), transform=ax.transAxes, ha='center', va='top')
ax.text(0.025, 0.95, 'b', transform=ax.transAxes, ha='left', va='top', fontweight='bold', fontsize=14)

# Save to file
plt.savefig('/home/atom/ongoing/work_stderr_dem/figures/final/Figure_S9_final.png', dpi=400)