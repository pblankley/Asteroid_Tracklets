# Imports
import os
import numpy as np
import scipy.interpolate
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import pandas as pd
import math
import kepcart as kc
import healpy as hp
import collections
import astropy
from collections import defaultdict
from collections import Counter
import MPC_library
import scipy.spatial
import pickle
from matplotlib.colors import Normalize

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

Observatories = MPC_library.Observatories

ObservatoryXYZ = Observatories.ObservatoryXYZ


## TODO: update docs and LOTS of paths for lower level functions below.  That information should
# be discussed with matt about file structure in our cohesive file and our demo file.

################## METHODS ##################################################
def is_two_line(line):
    """ This routine checks the 80-character input line to see if it contains a
    special character (S, R, or V) that indicates a 2-line record.
    ------
    Args: input line
    ------
    Returns: bool
    """
    note2 = line[14]
    return note2=='S' or note2=='R' or note2=='V'


def split_MPC_file(filename):
    """ This routine opens and reads filename, separating the records into those
    in the 1-line and 2-line formats. The 2-line format lines are merged into
    single 160-character records for processing line-by-line.
    -------
    Args: string; file path
    -------
    Returns: None; writes files to the directory you specify in filename argument
    """
    filename_1_line = filename.rstrip('.txt')+"_1_line.txt"
    filename_2_line = filename.rstrip('.txt')+"_2_line.txt"
    with open(filename_1_line, 'w') as f1_out, open(filename_2_line, 'w') as f2_out:
        line1=None
        with open(filename, 'r') as f:
            for line in f:
                if is_two_line(line):
                    line1=line
                    continue
                if line1 != None:
                    merged_lines = line1.rstrip('\n') + line
                    f2_out.write(merged_lines)
                    line1 = None
                else:
                    f1_out.write(line)
                    line1 = None

def get_sorted_tracklets(itf_filename):
    """ Takes in the path to the *.mpc file that contains tracklet information.
    ------
    Args: string; file path and name
    ------
    Returns: (dict, dict, list). Where the first dict, "tracklets" has the
        tracklet id as keys and the observation lines as values, the second dict
        "tracklets_jd_dict" has tracklet id as keys and julian data and value
        and the final list is the tracklet id's sorted by time from oldest
        to most recent.
    """
    tracklets = defaultdict(list)
    tracklets_jd_dict = {}
    with open(itf_filename) as infile:
        for line in infile:
            if not line.startswith('#'):
                desig = line[0:12]
                jd_tdb = float(line[43:57])
                if desig not in tracklets_jd_dict:
                    # Got a new tracklet
                    tracklets_jd_dict[desig] = jd_tdb
                tracklets[desig].append(line)
    sortedTrackletKeys = sorted(tracklets.keys(), key=lambda k: tracklets_jd_dict[k])
    return tracklets, tracklets_jd_dict, sortedTrackletKeys



def lunation_center(n, tref=2457722.0125, p=29.53055):
    """ Returns the jd of new moon, to the nearest half day. """
    t = tref + p*n
    tp = np.floor(t) + 0.5
    return tp

def separate_time_windows(tracklets, sortedTracklets, tracklets_jd_dict, \
                          file_stem,\
                          n_begin=-825, n_end=14, dt=15.):
    """ Sweep through the tracklets once, outputting them into a sequence of
    overlapping time ranges that can be processed separately. """
    t_center = lunation_center(n_begin)
    files = {}

    header='#trackletID yr   mn dy      obsCode mag filter  jd_tdb       x_target     y_target     z_target      x_obs       y_obs        z_obs     '

    for desig in sortedTracklets:
        jd_tdb = tracklets_jd_dict[desig]
        while(jd_tdb>t_center+dt):
            if n_begin in files:
                files[n_begin].close()
            n_begin +=1
            t_center = lunation_center(n_begin)
        for n in range(n_begin, n_end):
            if jd_tdb<lunation_center(n)-dt:
                break
            if n not in files:
                outfile = file_stem.rstrip('.mpc')+'_'+str(lunation_center(n))+'_pm'+str(dt)+'.mpc'
                files[n] = open(outfile, 'w')
                files[n].write(header+'\n')
            for line in tracklets[desig]:
                files[n].write(line)



def adjust_position(r, rho_target, re):
    """ This returns the topocentric distances and new heliocentric position
    vectors to the target, given the assumed distance "r" and the position
    vector of the observatory, "re".  This function transforms the data so the
    "view" of the asteroid is from the sun.
    -------
    Args: r; float, the assumed radial distance from the sun to the asteroid
          rho_target; a array with x,y,z vector from the mpc file with the
                    specified lunation center.
          re; a array with the x,y,z position vector for the observatory.
    -------
    Returns: (tuple, tuple) where the tuples return topocentric distances and
            the realted new heliocentric position vectors.
    """
    rho_x, rho_y, rho_z = rho_target
    xe, ye, ze = re
    Robs = np.sqrt(xe * xe + ye * ye + ze * ze)
    cos_phi = -(rho_x * xe + rho_y * ye + rho_z * ze) / Robs
    phi = np.arccos(cos_phi)
    sin_phi = np.sin(phi)

    xx2 = r*r - Robs*sin_phi * Robs*sin_phi

    if xx2 < 0:
        None, None

    xx = np.sqrt(xx2)
    yy = Robs * cos_phi

    rho_p = yy + xx

    # This could be done with numpy arrays
    x_p = xe + rho_p*rho_x
    y_p = ye + rho_p*rho_y
    z_p = ze + rho_p*rho_z

    rho_m = yy - xx

    # This could be done with numpy arrays
    x_m = xe + rho_m*rho_x
    y_m = ye + rho_m*rho_y
    z_m = ze + rho_m*rho_z

    return (rho_p, (x_p, y_p, z_p)), (rho_m, (x_m, y_m, z_m))



def transform_positions(n, r_func, file_stem, dt=45.):
    """ Does the transformations on the data using the date of the n-th new
    moon as the reference time. It is reading and processing the entire
    *.mpc file. This does the heliocentric tranformation for the assumed
    radius function, r_func.

    It then does light-time correction. And it appends a healpix number on
    each line in order to be able to quickly select data from a given
    region of sky.

    This generates a file called *.trans, and it incorporates the distance
    assumed in the file name.
    """
    infilename = file_stem.rstrip('.mpc')+'_'+str(lunation_center(n))+'_pm'+str(dt)+'.mpc'
    try:
      open(infilename, 'r')
    except IOError:
      return 0
    t_ref = lunation_center(n)
    r_ref = r_func(t_ref)
    r_name = "_r%.1lf" % (r_ref)
    outfilename = file_stem.rstrip('.mpc')+'_'+str(lunation_center(n))+'_pm'+str(dt)+r_name+'.trans'

    nside = 32

    #rot_mat = MPC_library.rotate_matrix(-MPC_library.Constants.ecl)

    with open(infilename, 'r') as infile, open(outfilename, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                header = line.rstrip()
                outfile.write(header + '          dt         x_cor       y_cor        z_cor       pix \n')
            else:
                lineID = line[:43]

                jd_tdb = float(line[43:57])

                x_target, y_target, z_target = line[57:97].split()
                r_target = np.array([float(x_target), float(y_target), float(z_target)])

                x_obs, y_obs, z_obs = line[97:135].split()
                r_obs = np.array([float(x_obs), float(y_obs), float(z_obs)])

                # Adjust positions
                dt = 0.0
                r_prev = r_func(jd_tdb-dt)
                rho_r_p, rho_r_m = adjust_position(r_prev, r_target, r_obs)
                dt = rho_r_p[0]/MPC_library.Constants.speed_of_light

                # Do light-time iterations.
                # Probably don't need to do this at this point, because it is
                # being re-done in a later step.
                i=0
                while(np.abs(r_func(jd_tdb-dt)-r_prev)>1e-8):
                    rho_r_p, rho_r_m = adjust_position(r_prev, r_target, r_obs)
                    dt = rho_r_p[0]/MPC_library.Constants.speed_of_light
                    r_prev = r_func(jd_tdb-dt)
                    i += 1

                x, y, z = rho_r_p[1]

                # At this point, I think I should fit a line
                # to the x, y, z components of each tracklet,
                # vs jd_tdb-dt-t_ref
                # rather than waiting for a later processing stage.
                # Then I'll have both a position and velocity vector.

                xp, yp, zp = x, y, z

                pix = hp.vec2pix(nside, xp, yp, zp, nest=True)

                outstring = line.rstrip() + " %13.6lf %12.7lf %12.7lf %12.7lf %5d\n"% \
                      (dt, xp, yp, zp, pix)

                outfile.write(outstring)


def xyz_to_proj_matrix(r_ref):
    """ This routine returns the 3-D rotation matrix for the given
    reference vector. """
    x_ref, y_ref, z_ref = r_ref
    r = np.sqrt(x_ref*x_ref + y_ref*y_ref + z_ref*z_ref)
    lon0 = np.arctan2(y_ref, x_ref)
    lat0 = np.arcsin(z_ref/r)
    slon0 = np.sin(lon0)
    clon0 = np.cos(lon0)
    slat0 = np.sin(lat0)
    clat0 = np.cos(lat0)

    mat = np.array([[-slon0, clon0, 0],
                    [-clon0*slat0, -slon0*slat0, clat0],
                    [clon0*clat0, slon0*clat0, slat0 ]])
    return mat

def select_positions_z(t_ref, z_ref, zdot_ref, vec, lines, outfilename):
    """
    # This is the one to use.  This routine will be used repeatedly.
    #
    # Trying a slightly different approach.
    # The set of lines that are being passed in have
    # been selected to be in the same region of sky
    # for an assumed distance.  We are going to re-transform
    # those assuming a fixed z (or gamma) value with respect
    # to the sun and the reference direction, rather than a
    # fixed r, at the reference time
    #
    # Rotate observatory positions to projection coordinates,
    # and recalculate simple z-based light-time correction.
    #
    # Rotate the observations to projection coordinates,
    # but they will be theta_x, theta_y only
    #
    # Fit the simple abg model, for fixed gamma and
    # possibly gamma_dot.
    """
    z_name = "_z%.1lf" % (z_ref)

    g = 1.0/z_ref
    gdot = zdot_ref/z_ref

    GM = MPC_library.Constants.GMsun

    # This rotation is taking things from equatorial to ecliptic
    rot_mat = MPC_library.rotate_matrix(-MPC_library.Constants.ecl)

    results_dict = defaultdict(list)

    # vec is the reference direction in equatorial coordinates
    # so we rotate to ecliptic, because we want to.
    vec = np.array(vec)
    ecl_vec = np.dot(rot_mat, vec)
    vec = ecl_vec
    vec = vec/np.linalg.norm(vec)
    # mat is a rotation matrix that converts from ecliptic
    # vectors to the projection coordinate system.
    # The projection coordinate system has z outward,
    # x parallel to increasing ecliptic longitude, and
    # y northward, making a right-handed system.
    mat = xyz_to_proj_matrix(vec)

    # Loop over all the lines from a *.trans file.
    for line in lines:
        if line.startswith('#'):
            # Concatenate all header lines?
            header = line.rstrip()
        else:
            lineID = line[:43]
            trackletID = line[0:12]

            jd_tdb = float(line[43:57])
            dtp = float(line[139:150])

            # Get unit vector to target
            x_target, y_target, z_target = line[58:98].split()
            r_target = np.array([float(x_target), float(y_target), float(z_target)])

            # Rotate to ecliptic coordinates
            r_target_ec = np.dot(rot_mat, r_target.T).T

            # Rotate to projection coordinates
            # Help me, Jesus!  Please use only one orientation
            # for vectors and matrices.
            theta_x, theta_y, theta_z = np.dot(mat, r_target_ec)

            # Ignore theta_z after this; it should be very nearly 1.

            # Get observatory position, ultimately in projection coordinates.
            x_obs, y_obs, z_obs = line[98:138].split()
            r_obs = np.array([float(x_obs), float(y_obs), float(z_obs)])

            # Rotate to ecliptic coordinates
            r_obs_ec = np.dot(rot_mat, r_obs.T).T

            # Rotate to projection coordinates
            xe, ye, ze = np.dot(mat, r_obs_ec)

            # This is the light travel time
            dlt = ze/MPC_library.Constants.speed_of_light

            # Append the resulting data to a dictionary keye do trackletID.
            results_dict[trackletID].append((jd_tdb, dlt, theta_x, theta_y, theta_z, xe, ye, ze))

    # Now that we have the observations for each tracklet gathered together,
    # we iterate through the tracklets, doing a fit for each one.
    results = []
    for k, v in results_dict.items():

        # Here's a version that incorporates radial gravitational
        # acceleration

        t_emit = [(obs[0]-obs[1]-t_ref) for obs in v]
        acc_z = -GM*g*g
        fac =[(1.0 + gdot*t + 0.5*g*acc_z*t*t - g*obs[7]) for obs, t in zip(v, t_emit)]

        A = np.vstack([t_emit, np.ones(len(t_emit))]).T

        # Can I put a simple acc_x term in here?
        x = [obs[2]*f + obs[5]*g for obs, f in zip(v, fac)]
        mx, cx = np.linalg.lstsq(A, x)[0]

        y = [obs[3]*f + obs[6]*g for obs, f in zip(v, fac)]
        my, cy = np.linalg.lstsq(A, y)[0]

        outstring = "%12s %16.9lf %16.9lf %16.9lf %16.9lf %16.9lf\n" % (k, cx, mx, cy, my, t_emit[0])
        results.append(outstring)

    if len(results)>0:
        with open(outfilename, 'w') as outfile:
            outstring = '#  z0 = %lf\n' % (z_ref)
            outfile.write(outstring)
            outstring = '#  zdot0 = %lf\n' % (zdot_ref)
            outfile.write(outstring)
            outstring = '#  vec= %lf, %lf, %lf\n' % (vec[0], vec[1], vec[2])
            outfile.write(outstring)
            outstring = '#  desig              alpha         alpha_dot       beta             beta_dot         dt \n'
            outfile.write(outstring)
            for outstring in results:
                outfile.write(outstring)



def select_positions_z_v2(t_ref, z_ref, zdot_ref, vec, lines, outfilename):
    """ # This is the one to use.  This routine will be used repeatedly.
    #
    # Trying a slightly different approach.
    # The set of lines that are being passed in have
    # been selected to be in the same region of sky
    # for an assumed distance.  We are going to re-transform
    # those assuming a fixed z (or gamma) value with respect
    # to the sun and the reference direction, rather than a
    # fixed r, at the reference time
    #
    # Rotate observatory positions to projection coordinates,
    # and recalculate simple z-based light-time correction.
    #
    # Rotate the observations to projection coordinates,
    # but they will be theta_x, theta_y only
    #
    # Fit the simple abg model, for fixed gamma and
    # possibly gamma_dot.
    """
    z_name = "_z%.1lf" % (z_ref)

    g = 1.0/z_ref
    gdot = zdot_ref/z_ref

    GM = MPC_library.Constants.GMsun

    rot_mat = MPC_library.rotate_matrix(-MPC_library.Constants.ecl)

    results_dict = defaultdict(list)

    vec = np.array(vec)
    ecl_vec = np.dot(rot_mat, vec)
    vec = ecl_vec
    vec = vec/np.linalg.norm(vec)
    mat = xyz_to_proj_matrix(vec)

    for line in lines:
        if line.startswith('#'):
            header = line.rstrip()
        else:
            lineID = line[:43]
            trackletID = line[0:12]

            jd_tdb = float(line[43:57])
            dtp = float(line[139:150])

            # Get unit vector to target
            x_target, y_target, z_target = line[58:98].split()
            r_target = np.array([float(x_target), float(y_target), float(z_target)])

            # Rotate to ecliptic coordinates
            r_target_ec = np.dot(rot_mat, r_target)

            # Rotate to projection coordinates
            theta_x, theta_y, theta_z = np.dot(mat, r_target_ec)

            # Ignore theta_z after this; it should be very nearly 1.

            # Get observatory position
            x_obs, y_obs, z_obs = line[98:138].split()
            r_obs = np.array([float(x_obs), float(y_obs), float(z_obs)])

            # Rotate to ecliptic coordinates
            r_obs_ec = np.dot(rot_mat, r_obs)

            # Rotate to projection coordinates
            xe, ye, ze = np.dot(mat, r_obs_ec)

            dlt = ze/MPC_library.Constants.speed_of_light

            results_dict[trackletID].append((jd_tdb, dlt, theta_x, theta_y, theta_z, xe, ye, ze))

    results = []
    for k, v in results_dict.items():

        # Here's a version that incorporates radial gravitational
        # acceleration

        # Solve this iteratively in this one loop section

        t_emit = [(obs[0]-obs[1]-t_ref) for obs in v]

        # We will approximate g_x(t), g_y(t), and g_z(t)
        # using a Taylor expansion.
        # The first two terms are zero by design.
        #
        # Given alpha, beta, gamma,
        # we would have r_0^2 = (alpha^2 + beta^2 + 1)*z_0^2
        # r^2 = (alpha^2 + beta^2 + 1)/gamma^2 ~ 1/gamma^2
        # g_x(t) ~ -0.5*GM*x_0*t^2/r^3
        # g_y(t) ~ -0.5*GM*y_0*t^2/r^3
        # g_z(t) ~ -0.5*GM*z_0*t^2/r^3
        #
        # We do not have alpha and beta initially,
        # but we assert gamma.
        #
        # We set alpha=beta=0 initially, least squares solve
        # the tracklets and obtain alpha, alpha-dot, beta,
        # and beta-dot.
        #
        # Then we use those values to estimate g_x,
        # g_y, and g_z for the next iteration.
        #
        # The process converges when alpha, alpha-dot,
        # beta, beta-dot do not change significantly.
        #
        # We could also do this same process with a
        # Kepler-stepper or a full n-body integration.

        alpha = beta = 0.0
        alpha_dot = beta_dot = 0.0
        cx, cy = 1.0, 1.0
        mx, my = 0.0, 0.0

        while(((cx-alpha)*(cx-alpha) + (cy-beta)*(cy-beta))>1e-16):

            alpha, beta = cx, cy
            alpha_dot, beta_dot = mx, my

            r2 = (alpha*alpha + beta*beta + 1.0)/(g*g)
            r3 = r2*np.sqrt(r2)
            r5 = r2*r3

            x0 = alpha/g
            y0 = beta/g
            z0 = 1.0/g

            vx0 = alpha_dot/g
            vy0 = beta_dot/g
            vz0 = gdot/g

            rrdot = x0*vx0 + y0*vy0 + z0*vz0

            acc_x = -GM*x0/r3
            acc_y = -GM*y0/r3
            acc_z = -GM*z0/r3

            jerk_x = -GM/r5*(r2*vx0 - 3.0*rrdot*x0)
            jerk_y = -GM/r5*(r2*vy0 - 3.0*rrdot*y0)
            jerk_z = -GM/r5*(r2*vz0 - 3.0*rrdot*z0)

            fac =[(1.0 + gdot*t + 0.5*g*acc_z*t*t + (1./6.0)*jerk_z*t*t*t - g*obs[7]) for obs, t in zip(v, t_emit)]

            A = np.vstack([t_emit, np.ones(len(t_emit))]).T

            x = [obs[2]*f + obs[5]*g - 0.5*g*acc_x*t*t - (1./6.0)*jerk_x*t*t*t for obs, f, t in zip(v, fac, t_emit)]
            mx, cx = np.linalg.lstsq(A, x)[0]

            y = [obs[3]*f + obs[6]*g - 0.5*g*acc_y*t*t - (1./6.0)*jerk_y*t*t*t for obs, f, t in zip(v, fac, t_emit)]
            my, cy = np.linalg.lstsq(A, y)[0]

        outstring = "%12s %16.9lf %16.9lf %16.9lf %16.9lf %16.9lf\n" % (k, cx, mx, cy, my, t_emit[0])
        results.append(outstring)

    if len(results)>0:
        with open(outfilename, 'w') as outfile:
            outstring = '#  z0 = %lf\n' % (z_ref)
            outfile.write(outstring)
            outstring = '#  zdot0 = %lf\n' % (zdot_ref)
            outfile.write(outstring)
            outstring = '#  vec= %lf, %lf, %lf\n' % (vec[0], vec[1], vec[2])
            outfile.write(outstring)
            outstring = '#  desig              alpha         alpha_dot       beta             beta_dot         dt \n'
            outfile.write(outstring)
            for outstring in results:
                outfile.write(outstring)




def cluster_positions_z(t_ref, z_zdot_pairs, vec, lines):
    """ # Here I am doing the same thing as the previous routine, but without files.
    #
    # It takes a reference time (t_ref), a set of z, zdot pairs (z_zdot_pairs),
    # a reference direction vector (vec), and a set of observation lines that
    # have been selected for a region of sky and time slice (lines)
    #
    # It returns a dictionary of results that have z, zdot pairs as keys and
    # sets of fitted tracklets as results.  Each result has the form:
    #
    # trackletID alpha alpha_dot beta beta_dot t_emit,
    # where t_emit is the light time-corrected time relative to the reference
    # time.  The coordinates are now in tangent plane projection.
    """
    GM = MPC_library.Constants.GMsun

    rot_mat = MPC_library.rotate_matrix(-MPC_library.Constants.ecl)

    results_dict = defaultdict(list)

    vec = np.array(vec)
    ecl_vec = np.dot(rot_mat, vec.T).T
    vec = ecl_vec
    vec = vec/np.linalg.norm(vec)
    mat = xyz_to_proj_matrix(vec)

    for line in lines:
        if line.startswith('#'):
            header = line.rstrip()
        else:
            lineID = line[:43]
            trackletID = line[0:12]

            jd_tdb = float(line[43:57])
            dtp = float(line[139:150])

            # Get unit vector to target
            x_target, y_target, z_target = line[58:98].split()
            r_target = np.array([float(x_target), float(y_target), float(z_target)])

            # Rotate to ecliptic coordinates
            r_target_ec = np.dot(rot_mat, r_target)

            # Rotate to projection coordinates
            theta_x, theta_y, theta_z = np.dot(mat, r_target_ec)

            # Ignore theta_z after this; it should be very nearly 1.

            # Get observatory position
            x_obs, y_obs, z_obs = line[98:138].split()
            r_obs = np.array([float(x_obs), float(y_obs), float(z_obs)])

            # Rotate to ecliptic coordinates
            r_obs_ec = np.dot(rot_mat, r_obs.T).T

            # Rotate to projection coordinates
            xe, ye, ze = np.dot(mat, r_obs_ec)

            dlt = ze/MPC_library.Constants.speed_of_light

            results_dict[trackletID].append((jd_tdb, dlt, theta_x, theta_y, theta_z, xe, ye, ze))

    # All the work done above is independent of the z0 and zdot0 values

    master_results = {}
    for z_zdot in z_zdot_pairs:
        #print(z_zdot)
        z, zdot = z_zdot
        g = 1.0/z
        gdot = zdot/z
        results = []
        for k, v in results_dict.items():

            # Here's a version that incorporates radial gravitational
            # acceleration

            t_emit = [(obs[0]-obs[1]-t_ref) for obs in v]
            acc_z = -GM*g*g
            fac =[(1.0 + gdot*t + 0.5*g*acc_z*t*t - g*obs[7]) for obs, t in zip(v, t_emit)]

            A = np.vstack([t_emit, np.ones(len(t_emit))]).T

            x = [obs[2]*f + obs[5]*g for obs, f in zip(v, fac)]
            mx, cx = np.linalg.lstsq(A, x)[0]

            y = [obs[3]*f + obs[6]*g for obs, f in zip(v, fac)]
            my, cy = np.linalg.lstsq(A, y)[0]

            result = (k, cx, mx, cy, my, t_emit[0])
            results.append(result)

        master_results[z_zdot] = results

    return master_results

def cluster_sky_regions(z_zdot_pairs, pixels, infilename='/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line_2457397.5_pm15.0_r2.5.trans', nside=8, n=-11, angDeg=5.5):
    hp_dict = defaultdict(list)
    with open(infilename) as file:
        for line in file:
            if line.startswith('#'):
                continue
            pix = int(line.split()[-1])
            hp_dict[pix].append(line)

    pixel_results = {}
    #for i in range(hp.nside2npix(nside)):
    for i in pixels:
        print(i)
        # Probably don't need to repeat the vector neighbor calculation.
        # This can just be stored.
        vec = hp.pix2vec(nside, i, nest=True)
        neighbors = hp.query_disc(32, vec, angDeg*np.pi/180., inclusive=True, nest=True)
        lines = []
        for pix in neighbors:
            for line in hp_dict[pix]:
                lines.append(line)
        if len(lines) > 0:
            #print(i, len(lines))
            pixel_results[i] = cluster_positions_z(lunation_center(n), z_zdot_pairs, vec, lines)

    return pixel_results

def member_counts(k, sep='|', suff='_'):
    keys = k.split(sep)
    stems = [key.split(suff)[0] for key in keys]
    stem_counter = Counter(stems)
    return stem_counter

def unique_clusters(test_set):
    check_dict = {}
    errs=[]
    for k in test_set:
        stem_counter = member_counts(k)
        if len(stem_counter)>1:
            errs.append(k)
        else:
            for k, v in stem_counter.items():
                if k not in check_dict:
                    check_dict[k] = v
                elif v > check_dict[k]:
                    check_dict[k] = v
    return check_dict, errs


zs = np.arange(2.5, 2.6, 0.5)
#zs = np.arange(3.33, 3.5, 0.5)
#zs = np.arange(10.0, 10.1, 0.5)
zdots = np.arange(-1e-2, 1.1e-2, 2.0e-3)
z_zdots = [(x,y) for x in zs for y in zdots]

def do_training_run(pixels, z_zdots=z_zdots, infilename='/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_2457397.5_pm15.0_r2.5.trans'):
    master = cluster_sky_regions(z_zdots, pixels, infilename=infilename)

    dt_range = np.arange(10, 80, 10)
    results_dict = {}
    # In the process of clustering here, we need to
    # keep track of the 'best' or biggest cluster that
    # each tracklet is a member of.
    rates_dict={}
    for dt in dt_range:
        for rad in np.arange(0.0001, 0.0100, 0.0001):
            cluster_counter = Counter()
            for pix, d in master.items():
                for z_zdot, arrows in d.items():
                    i = 0
                    label_dict={}
                    combined=[]
                    for k, cx, mx, cy, my, t in arrows:
                        label_dict[i] = k
                        combined.append([cx, mx*dt, cy, my*dt])
                        i +=1
                    points=np.array(combined)
                    tree = scipy.spatial.cKDTree(points)
                    matches = tree.query_ball_tree(tree, rad)
                    for j, match in enumerate(matches):
                        if len(match)>1:
                            cluster_list =[]
                            tracklet_params=[]
                            for idx in match:
                                cluster_list.append(label_dict[idx].strip())
                                tracklet_params.append(combined[idx])
                            cluster_key='|'.join(sorted(cluster_list))
                            cluster_counter.update({cluster_key: 1})

            errs = 0
            for i, k in enumerate(cluster_counter.keys()):
                keys = k.split('|')
                stems = [key.split('_')[0] for key in keys]
                stem_counter = Counter(stems)
                if len(stem_counter)>1:
                    errs +=1

            rates_dict[dt, rad] = cluster_counter.keys(), errs

    for dt in dt_range:
        values = []
        for k, v in rates_dict.items():
            dtp, d = k
            if dtp==dt:
                test_set = list(v[0])
                ncs, nes = len(unique_clusters(test_set)[0]), len(unique_clusters(test_set)[1])
                values.append((d, ncs, nes))

        values = sorted(values, key=lambda v: v[0])
        ds = [v[0] for v in values]
        nclusters = [v[1] for v in values]
        nerrors = [v[2] for v in values]
        results_dict[dt] = ds, nclusters, nerrors

    return results_dict


def output_sky_regions(z_zdot_pairs, pixels, infilename='/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line_2457397.5_pm45.0_r2.5.trans', nside=8, n=-11, angDeg=5.5):
    hp_dict = defaultdict(list)
    with open(infilename) as file:
        i=0
        for line in file:
            if line.startswith('#'):
                continue
            pix = int(line.split()[-1])
            hp_dict[pix].append(line)
            i += 1

    pixel_results = {}
    #for i in range(hp.nside2npix(nside)):
    for i in pixels:
        # Probably don't need to repeat the vector neighbor calculation.
        # This can just be stored.
        vec = hp.pix2vec(nside, i, nest=True)
        neighbors = hp.query_disc(32, vec, angDeg*np.pi/180., inclusive=True, nest=True)
        lines = []
        for pix in neighbors:
            for line in hp_dict[pix]:
                lines.append(line)

    return lines



def accessible_clusters(pixels, infilename='/Users/paulblankley/Desktop/CS 182/Project/data/UnnObs_Training_1_line_A_2457397.5_pm15.0_r2.5.trans', z_zdots=z_zdots):
    true_counts={}
    for pix in pixels:
        lines = output_sky_regions(z_zdots, [pix], infilename=infilename)
        #print(len(lines))
        trackletCounter = Counter()
        for line in lines:
            trackletID=line.split()[0]
            trackletCounter.update({trackletID : 1})

        mergedCounter = Counter()
        for k, v in trackletCounter.items():
            mergedCounter.update({k[:-4]:1})
        true_counts[pix]=len([k for k, v in mergedCounter.items() if v>1])

    return true_counts

def do_run(z_zdot_pairs, pixels, infilename, nside=8, n=-11, angDeg=5.5, dt=50, rad=0.0008):

    print('tranforming')
    master=cluster_sky_regions(z_zdot_pairs, pixels, infilename=infilename, nside=nside, n=n, angDeg=angDeg)

    print('clustering')
    #results_dict = {}
    # In the process of clustering here, we need to
    # keep track of the 'best' or biggest cluster that
    # each tracklet is a member of.
    #rates_dict={}
    cluster_counter = Counter()
    cluster_dict={}
    for pix, d in master.items():
        for z_zdot, arrows in d.items():
            i = 0
            label_dict={}
            combined=[]
            for k, cx, mx, cy, my, t in arrows:
                label_dict[i] = k
                combined.append([cx, mx*dt, cy, my*dt])
                i +=1
            points=np.array(combined)
            tree = scipy.spatial.cKDTree(points)
            matches = tree.query_ball_tree(tree, rad)
            for j, match in enumerate(matches):
                if len(match)>1:
                    cluster_list =[]
                    tracklet_params=[]
                    for idx in match:
                        cluster_list.append(label_dict[idx].strip())
                        tracklet_params.append(combined[idx])
                    cluster_key='|'.join(sorted(cluster_list))
                    cluster_counter.update({cluster_key: 1})

    return cluster_counter, cluster_dict

def get_original_tracklets_dict(itf_filename='/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line.txt'):
    tracklets = defaultdict(list)
    with open(itf_filename) as infile:
        for line in infile:
            if not line.startswith('#'):
                desig = line[0:12].strip()
                tracklets[desig].append(line)
    return tracklets

def get_original_data(cluster_key, tracklets_dict):
    lines=[]
    trackletIDs = cluster_key.split('|')
    for trackletID in trackletIDs:
        for line in tracklets_dict[trackletID]:
            lines.append(line)
    lines = sorted(lines, key=lambda x: x[15:31])
    return lines

def generate_sky_region_files(infilename='/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line_2457397.5_pm15.0_r2.5.trans', nside=8, n=-11, angDeg=5.5, z0=2.5, zdot0=0.0):
    hp_dict = defaultdict(list)
    with open(infilename) as file:
        for line in file:
            if line.startswith('#'):
                continue
            pix = int(line.split()[-1])
            hp_dict[pix].append(line)

    for i in range(hp.nside2npix(nside)):
        vec = hp.pix2vec(nside, i, nest=True)
        neighbors = hp.query_disc(32, vec, angDeg*np.pi/180., inclusive=True, nest=True)
        lines = []
        for pix in neighbors:
            for line in hp_dict[pix]:
                lines.append(line)
        outfilename = infilename.rstrip('.trans') + '_hp_' + ('%03d' % (i)) + '_z'+ ('%.2lf' % (z0))+'_zdot' + ('%+5.1le' % (zdot0))
        if len(lines) > 0:
            select_positions_z(lunation_center(n), z0, zdot0, vec, lines, outfilename)

def generate_sky_region_files_v2(infilename='/Users/paulblankley/Desktop/CS 182/Project/data/itf_new_1_line_2457397.5_pm15.0_r2.5.trans', nside=8, n=-11, angDeg=5.5, z0=2.5, zdot0=0.0):
    hp_dict = defaultdict(list)
    with open(infilename) as file:
        for line in file:
            if line.startswith('#'):
                continue
            pix = int(line.split()[-1])
            hp_dict[pix].append(line)

    for i in range(hp.nside2npix(nside)):
        vec = hp.pix2vec(nside, i, nest=True)
        neighbors = hp.query_disc(32, vec, angDeg*np.pi/180., inclusive=True, nest=True)
        lines = []
        for pix in neighbors:
            for line in hp_dict[pix]:
                lines.append(line)
        outfilename = infilename.rstrip('.trans') + '_hp_' + ('%03d' % (i)) + '_z'+ ('%.2lf' % (z0))+'_zdot' + ('%+5.1le_v2' % (zdot0))
        if len(lines) > 0:
            select_positions_z_v2(lunation_center(n), z0, zdot0, vec, lines, outfilename)



def make_figure(filename):
    plt.ioff()
    mxs, cxs, mys, cys, dts =[], [], [], [], []
    for line in open(filename):
        if line.startswith('#'):
            continue
        desig, cx, mx, cy, my, dt = line.split()
        mxs.append(float(mx))
        cxs.append(float(cx))
        mys.append(float(my))
        cys.append(float(cy))
        dts.append(float(dt))

    fig=plt.figure(figsize=(18, 16))

    #norm = Normalize()
    #norm.autoscale(colors)

    colormap = cm.inferno

    plt.quiver(cxs, cys, mxs, mys, dts, scale=0.3, width=0.0003)

    plt.xlim(-0.2, 0.2)
    plt.ylim(-0.2, 0.2)
    plt.xlabel('alpha')
    plt.ylabel('beta')
    outfile = filename+'.pdf'
    plt.savefig(outfile)
    plt.close()
    plt.ion()
