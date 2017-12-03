# Imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import utils  # Get the module structured so this works.

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
# os.path.join(BASE_DIR, 'data/plist_df.json')

## TODO: add __init__.py file to the itf directory and make the whole package more module-like
## TODO: ask matt about orbit fitting process and  about what we want the output of the code to be
## QUESTION: Do we want the output to be the percent of confirmed (by orbit fitting) correct culstering
            # matches accompanied by a cool colored graph, or are we looking to output something else?

""" Here we will have a demo for the code that gets a subset of data from the ITF
    and runs our clustering process over the subset. It will then give the
    resulting clusters and the realted orbit-fit validated results, along with
    the total number of asteroids we accurately matched, based on the orbit fitting.
"""


mpc_path = '/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line.mpc'
training_path = '/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A.mpc'

###############################################################################

tracklets, tracklets_jd_dict, sortedTracklets = get_sorted_tracklets(mpc_path)
UnnObs_tracklets, UnnObs_tracklets_jd_dict, UnnObs_sortedTracklets = get_sorted_tracklets(training_path)

# In the +-15 day range from the lunar cycle
separate_time_windows(tracklets, sortedTracklets, tracklets_jd_dict, file_stem=mpc_path, dt=15.)
separate_time_windows(UnnObs_tracklets, UnnObs_sortedTracklets, UnnObs_tracklets_jd_dict, file_stem=training_path, dt=15.)

# In the +-45 day range from the lunar cycle (only training data)
separate_time_windows(UnnObs_tracklets, UnnObs_sortedTracklets, UnnObs_tracklets_jd_dict, file_stem=training_path, dt=45.)

# Transform the positions once per dictance class for +-15 day range from the new moon
for n in range(-825,14):
    transform_positions(n, lambda t: 2.5, file_stem=mpc_path, dt=15.)

# Transform the positions once per distance class for +-15 and +-45 for training data
for n in range(-825,14):
    transform_positions(n, lambda t: 2.5, file_stem=training_path, dt=15.)
    transform_positions(n, lambda t: 2.5, file_stem=training_path, dt=45.)

zs = np.arange(2.5, 2.6, 0.5)
#zs = np.arange(3.33, 3.5, 0.5)
#zs = np.arange(10.0, 10.1, 0.5)
zdots = np.arange(-1e-2, 1.1e-2, 2.0e-3)
z_zdots = [(x,y) for x in zs for y in zdots]

pix_runs = {}
nside=8

for pix in range(hp.nside2npix(nside)):
    pix_runs[pix] = do_training_run([pix], z_zdots=z_zdots,
        infilename='/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_2457397.5_pm15.0_r2.5.trans')

with open('UnnObs_pix_runs_2457397.5_pm15.0_z2.5_v2.pickle', 'wb') as handle:
    pickle.dump(pix_runs, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open('UnnObs_pix_runs_2457397.5_pm15.0_z2.5_v2.pickle', 'rb') as handle:
    pix_runs_b = pickle.load(handle)

assert(pix_runs == pix_runs_b)

UnnObs378_lines = output_sky_regions(z_zdots, [378], infilename='/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_2457397.5_pm45.0_r2.5.trans')


true_count_dict=accessible_clusters(list(pix_runs.keys()))
true_count=sum(true_count_dict.values())

with open('UnnObs_pix_runs_2457397.5_pm15.0_z2.5_v2.pickle', 'rb') as handle:
    pix_runs = pickle.load(handle)

from operator import add

# Deciding on a distance value for cluster radius with different scales for
# position vs velocity importance
for dt in np.arange(10, 80, 10):
    pixels=list(pix_runs.keys())
    ds = pix_runs[pixels[0]][dt][0]
    nclusters = pix_runs[pixels[0]][dt][1]
    nerrors = pix_runs[pixels[0]][dt][2]
    for pix in pixels[1:]:
        nclusters = list(map(add, nclusters, pix_runs[pix][dt][1]))
        nerrors = list(map(add, nerrors, pix_runs[pix][dt][2]))
    nclusters=np.array(nclusters)

    plt.plot(ds, nclusters, label=dt)

plt.axhline(true_count, ls='dashed')
plt.xlabel('d (cluster radius)')
plt.ylabel('N clusters')
plt.text(0.005, 400, r'$\gamma=0.4$', fontsize=15)
plt.legend()
plt.savefig('nclusters_vs_d_z2.5_v2.pdf')
plt.show()

## number of errors for the same above metrics
for dt in np.arange(10, 80, 10):
    pixels=list(pix_runs.keys())
    ds = pix_runs[pixels[0]][dt][0]
    nclusters = pix_runs[pixels[0]][dt][1]
    nerrors = pix_runs[pixels[0]][dt][2]
    for pix in pixels[1:]:
        nclusters = list(map(add, nclusters, pix_runs[pix][dt][1]))
        nerrors = list(map(add, nerrors, pix_runs[pix][dt][2]))
    nclusters=np.array(nclusters)
    nerrors=np.array(nerrors)
    '''
    if dt==50:
        for d, nc, ne in zip(ds, nclusters, nerrors):
            print(d, nc, ne)
    '''


    plt.plot(ds, nerrors, label=dt)

plt.ylim((0,3000))
plt.xlabel('d (cluster radius)')
plt.ylabel('N errors')
plt.text(0.003, 2000, r'$\gamma=0.4$', fontsize=15)
plt.legend()
plt.savefig('nerrors_vs_d_z2.5.pdf')
plt.show()
#

# AUC Curves for the different velo/position weights
ntrue=sum(true_count_dict.values())
for dt in np.arange(10, 80, 10):
    pixels=list(pix_runs.keys())
    ds = pix_runs[pixels[0]][dt][0]
    nclusters = pix_runs[pixels[0]][dt][1]
    nerrors = pix_runs[pixels[0]][dt][2]
    for pix in pixels[1:]:
        nclusters = list(map(add, nclusters, pix_runs[pix][dt][1]))
        nerrors = list(map(add, nerrors, pix_runs[pix][dt][2]))
    nclusters=np.array(nclusters)
    nerrors=np.array(nerrors)

    plt.plot(nerrors/ntrue, nclusters/ntrue, label=dt)

plt.xlim((0,0.1))
plt.ylim((0, 1))
plt.xlabel('Error rate')
plt.ylabel('Fraction correct')
plt.text(0.05, 0.2, r'$\gamma=0.4$', fontsize=15)
plt.legend()
plt.savefig('AUCish_z2.5.pdf')
plt.show()
#

# Do the actual run
zs = np.arange(2.5, 2.6, 0.5)
zdots = np.arange(-1e-2, 1.1e-2, 2.0e-3)
z_zdots = [(x,y) for x in zs for y in zdots]
nside=8
pixels= range(hp.nside2npix(nside))

thepath = '/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line_2457397.5_pm15.0_r2.5.trans'
first_real_run = do_run(z_zdots, pixels, infilename=thepath, nside=8, n=-11, angDeg=5.5, dt=50, rad=0.0008)

with open('itf_new_1_line_2457397.5_pm15.0_z2.5.pickle', 'wb') as handle:
    pickle.dump(first_real_run, handle, protocol=pickle.HIGHEST_PROTOCOL)


with open('itf_new_1_line_2457397.5_pm15.0_z2.5.pickle', 'rb') as handle:
    first_real_run_b = pickle.load(handle)

assert(first_real_run == first_real_run_b)

original_tracklets_dict = get_original_tracklets_dict()

with open('obs.obs', 'w') as file:
    for cluster_key in first_keys:
        file.write("%s\n" %(cluster_key))
        lines=get_original_data(cluster_key, original_tracklets_dict)
        for line in lines:
            file.write(line)

generate_sky_region_files()

z0=3.0
for n in range(-11, -9):
    infilename = '/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line_'+str(lunation_center(n))+'_pm15.0_r2.5.trans'
    print(infilename)
    generate_sky_region_files(infilename=infilename, n=n, z0=z0)

zdot0=+2.5e-3
for z0 in [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5]:
    for n in range(-11, -10):
        infilename = '/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_'+str(lunation_center(n))+'_pm45.0_r2.5.trans'
        print(infilename)
        generate_sky_region_files_v2(infilename=infilename, n=n, z0=z0, zdot0=zdot0)

z0=3.1
#for zdot0 in [-1.0e-2, -8.0e-3, -6.0e-3, -4.0e-3, -2.0e-3, 0.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2]:
for zdot0 in [-4.0e-3]:
    for n in range(-11, -10):
        infilename = '/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_'+str(lunation_center(n))+'_pm45.0_r2.5.trans'
        print(infilename)
        generate_sky_region_files_v2(infilename=infilename, n=n, z0=z0, zdot0=zdot0)

z0=3.1
#for zdot0 in [-1.0e-2, -8.0e-3, -6.0e-3, -4.0e-3, -2.0e-3, 0.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3, 1.0e-2]:
for zdot0 in [-4.0e-3]:
    for n in range(-11, -10):
        infilename = '/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_'+str(lunation_center(n))+'_pm45.0_r2.5.trans'
    for n in range(-11, -10):
        infilename = '/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_'+str(lunation_center(n))+'_pm45.0_r2.5.trans'
        filename = infilename.rstrip('.trans') + '_hp_378' + ('_z%4.2lf_zdot%+3.1le_v2' % (z0, zdot0))
        make_figure(filename)


# As many 'make_figure' calls as needed
make_figure('/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_2457397.5_pm45.0_r2.5_hp_378_z2.00_zdot+2.5e-03_v2')
