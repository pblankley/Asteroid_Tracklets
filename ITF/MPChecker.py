#!/usr/bin/env python
# -*- coding: utf-8 -*-

# I need to set this up so that we call it
# with a file that has a set of exposures we want to
# check (all from the same day), the HP dictionary
# for that day, and the orbit file that covers that
# day.

# This was save from the MPC.ipynb notebook.  15 Oct 2017
# I am modifying it to work with new orbit and healpix dictionaries
#
## MPChecker
#### Matthew J. Holman

# 11 June 2016
# 
# This is supposed to be just the routines from my MPC notebook that are relevant to MPChecker.
# 
# You will need to make sure you have a copy of the kepcart library. 
#
# You will also need:
#
# The tai-utc.dat file
# http://maia.usno.navy.mil/ser7/tai-utc.dat
#
# The list of observatory codes at
# http://www.minorplanetcenter.net/iau/lists/ObsCodes.html
# saved as ObsCodes.txt
#
# For dealing with leap seconds and polar motion
# finals2000A.all at
# ftp://maia.usno.navy.mil/ser7/finals2000A.all

import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import math
import kepcart as kc
import healpy as hp
import collections
import pickle
import time
import scipy.interpolate
import sys
import MPC_library as MPC
import ele220
from os import listdir
from os.path import isfile, join
import h5py

class Constants:
    GMsun = 2.9591220828559115e-04 
    Rearth_km = 6378.1363
    au_km = 149597870.700 # This is now a definition 
    Rearth_AU = Rearth_km/au_km
    ecl = (84381.4118*(1./3600)*np.pi/180.) # Obliquity of ecliptic at J2000
    speed_of_light = 2.99792458e5 * 86400./au_km

# ### The NOVAS package
# 
# First, let's get the USNO's python NOVAS package.  We'll need that.
# 
# http://aa.usno.navy.mil/software/novas/novas_py/novaspy_intro.php
# 
# Just type 
# 
# pip install novas
# 
# pip install novas_de405
# 
# Here's the reference:
# 
# Barron, E. G., Kaplan, G. H., Bangert, J., Bartlett, J. L., Puatua, W., Harris, W., & Barrett, P. (2011)
# “Naval Observatory Vector Astrometry Software (NOVAS) Version 3.1, Introducing a Python Edition,” Bull. AAS, 43, 2011.

from novas import compat as novas
from novas.compat import eph_manager
from novas.compat import solsys
# This opens DE405 by default.
jd_start, jd_end, number = eph_manager.ephem_open()

class Observatory:

    # Parses a line from the MPC's ObsCode.txt file
    def parseObsCode(self, line):
        code, longitude, rhocos, rhosin, ObsName = line[0:3], line[4:13], line[13:21], line[21:30], line[30:].rstrip('\n')
        if longitude.isspace():
            longitude = None
        if rhocos.isspace():
            rhocos = None
        if rhosin.isspace():
            rhosin = None
        return code, longitude, rhocos, rhosin, ObsName

    def __init__(self):

        self.observatoryPositionCache = {} # previously calculated positions to speed up the process

        # Convert ObsCodes.txt lines to geocentric x,y,z positions and
        # store them in a dictionary.  The keys are the observatory
        # code strings, and the values are (x,y,z) tuples.
        # Spacecraft and other moving observatories have (None,None,None)
        # as position.
        ObservatoryXYZ = {}
        with open('ObsCodes.txt', 'r') as f:
            next(f)
            for line in f:
                code, longitude, rhocos, rhosin, Obsname = self.parseObsCode(line)
                if longitude and rhocos and rhosin:
                    rhocos, rhosin, longitude = float(rhocos), float(rhosin), float(longitude)
                    longitude *= np.pi/180.
                    x = rhocos*np.cos(longitude)
                    y = rhocos*np.sin(longitude)
                    z = rhosin
                    ObservatoryXYZ[code]=(x,y,z)
                else:
                    ObservatoryXYZ[code]=(None,None,None)
        self.ObservatoryXYZ = ObservatoryXYZ

    # The routine below calculates the heliocentric position of the observatory
    # in equatorial cartesian coordinates.
    def getObservatoryPosition(self, obsCode, jd_utc):

        # If obsCode==None, set the observatory to be the geocenter
        if not obsCode:
            obsCode = '500'

        if (obsCode, jd_utc) in self.observatoryPositionCache:
            return self.observatoryPositionCache[(obsCode, jd_utc)]
    
        obsVec = self.ObservatoryXYZ[obsCode]

        if obsVec[0]==None:
            print "problem with obscode: ", obsCode
            return None, None, None
        else:

            jd_tdb  = EOP.jdTDB(jd_utc)
            pos = getEarthPosition(jd_tdb)

            if obsCode=='500':
                geocentric_vec = np.zeros(3)
            else:
                delta_t = EOP.delta_t(jd_utc)
                jd_ut1  = EOP.jdUT1(jd_utc)
                xp, yp  = EOP.pmx(jd_utc), EOP.pmy(jd_utc)

                geocentric_vec = Constants.Rearth_AU*np.array(novas.ter2cel(jd_ut1, 0.0, delta_t, xp, yp, obsVec))

            heliocentric_vec = pos+geocentric_vec
            self.observatoryPositionCache[(obsCode, jd_utc)] = heliocentric_vec
            return heliocentric_vec

Observatories = Observatory()
ObservatoryXYZ = Observatories.ObservatoryXYZ


class EarthAndTime:
    # Dealing with leap seconds and polar motion
    # ### Relating this to MPC data
    # 
    # I believe that the MPC observations have dates in UTC, which is the conventional thing to do.   
    # According to Gareth Williams, prior to 1972 Jan 1 the times are probably UT1.
    # 
    # So we take a time from an MPC observation.  If it's prior to 1972 Jan 1 we assume it's UT1 and 
    # we get delta_t (TT-UT1) from the historical table.  If it's on or after 1972 Jan 1 we determine 
    # the number of leap seconds from function below and then calculate delta_t.  
    # 
    def __init__(self, filename1='finals2000A.all', filename2='tai-utc.dat'):
        _xydeltat = {}
        with open(filename1) as f:
            _mjd    = slice(7, 15)
            _UT1_UTC= slice(58, 68)
            _pmx    = slice(18, 27)
            _pmy    = slice(37, 46)
    
            for line in f:
                if not line[_UT1_UTC].strip() == '':
                    _jd     = float(line[_mjd]) + 2400000.5
                    UT1_UTC= float(line[_UT1_UTC])
                    pmx    = float(line[_pmx])
                    pmy    = float(line[_pmy])
                    _xydeltat[_jd] = (UT1_UTC, pmx, pmy)

        jds = sorted(_xydeltat.keys())
        ut1_utcs = [_xydeltat[jd][0] for jd in jds]
        pmxs     = [_xydeltat[jd][1] for jd in jds]
        pmys     = [_xydeltat[jd][2] for jd in jds]

        self.ut1_utc_func = scipy.interpolate.interp1d(jds, ut1_utcs)
        self.pmx_func     = scipy.interpolate.interp1d(jds, pmxs)
        self.pmy_func     = scipy.interpolate.interp1d(jds, pmys)

        # Get the TAI-UTC data from:
        # http://maia.usno.navy.mil/ser7/tai-utc.dat
        self.tai_minus_utc_dict = {}
        with open(filename2) as f:
            _jd      = slice(16, 28)
            _tai_minus_utc = slice(36, 49)
            _tref    = slice(59, 66)
            _coeff   = slice(69, 79)

            for line in f:
                tai_minus_utc = float(line[_tai_minus_utc])
                jd      = float(line[_jd])
                tref    = float(line[_tref])
                coeff   = float(line[_coeff])
                self.tai_minus_utc_dict[jd] = tai_minus_utc, tref, coeff

    def pmx(self, jd_utc):
        return self.pmx_func(jd_utc)

    def pmy(self, jd_utc):
        return self.pmy_func(jd_utc)

    # TDT = UTC + (TAI-UTC) + 32.184 sec
    # UT1 = UTC + delta_t
    # TT - TDB is less than 2 milliseconds.

    def ut1_utc(self, jd_utc):
        return self.ut1_utc_func(jd_utc)

    def tai_utc(self, jd_utc):
        if jd_utc<2437300.5:
            return jd_utc, 0.0
        ks = sorted(self.tai_minus_utc_dict.keys())
        m = max(i for i in ks if (i-jd_utc)<0.0)
        base, tref, coeff = self.tai_minus_utc_dict[m]
        tai_m_utc = base + coeff*(jd_utc-240000.5-tref)
        return tai_m_utc
    
    def jdTT(self, jd_utc):
        leaps = self.tai_utc(jd_utc)
        jd_tt = jd_utc + (32.184 + leaps)/(24.0*60*60)
        return jd_tt

    def jdUT1(self, jd_utc):
        DUT1  = self.ut1_utc(jd_utc)
        jd_ut1 = jd_utc + DUT1/(24.0*60*60)
        return jd_ut1

    def jdTDB(self, jd_utc):
        jd_tt = self.jdTT(jd_utc)
        _, tdb_tt = novas.tdb2tt(jd_tt)
        jd_tdb = jd_tt + tdb_tt/(24.0*60*60)
        return jd_tdb

    def delta_t(self, jd_utc):
        leaps = self.tai_utc(jd_utc)
        DUT1  = self.ut1_utc(jd_utc)
        delta_t = 32.184 + leaps - DUT1
        return delta_t

EOP = EarthAndTime()

def getEarthPosition(jd_tdb):
    pos, _ = solsys.solarsystem(jd_tdb, 3, 1)
    return pos

def parseOrbit(self, GMsun, line):
    (desig, H, G, q, e, incl, longnode, argperi, tperi, epoch, _, _, _) = line.split()
    elts = (q, e, incl, longnode, argperi, tperi, epoch)
    return H, G, elts


def readMPCorbits(filename='mpn_K12CJ.human'):
    with open(filename) as fp:
        orbit_dict = pickle.load(fp)
    return orbit_dict

def readMPCfile(filename='mpn_K12CJ.txt'):
    orb_dict = {}
    with open(filename) as fd:
        for i, line in enumerate(fd):
            e = ele220.Ele220(line)
            orb_dict[e.desig()] = e.h(), e.g(), (e.periDist(), e.ecc(), e.inc()*np.pi/180., e.node()*np.pi/180., e.argPeri()*np.pi/180., e.timePeri())
    return orb_dict

def readMPC_HDF5_file(filename='mpn_K12CJ.human.hdf5'):
    orb_dict = {}
    with h5py.File(filename, 'r') as fd:
        data = fd['orbit'][()]
        print data.shape, type(data), type(data[0])
        keys = data[:, 0]
        values = data[:, 1:]
        orb_dict = dict(zip(keys, values))
        '''
        for row in data:
            desig, H, G, q, e, incl, longnode, argperi, tperi, epoch, u, rms = row
            print desig, rms
            orb_dict[desig] = H, G, (q, e, incl*np.pi/180., longnode*np.pi/180., argperi*np.pi/180., tperi)
        '''
    return orb_dict

def readPixels(filename):
    with open(filename) as fp:
        pixel_dict = pickle.load(fp)
    return pixel_dict

def readHEALPix(filename):
    hp_dict = collections.defaultdict(list)
    with open(filename) as fp:
        for line in fp:
            desig, ra, dec, mag, hp = line.split()
            hp = int(hp)
            hp_dict[hp].append(desig)
    return hp_dict

def getPVs(GMsun, orbits, jd_tt, lt):
    q        = orbits[:,0].copy()
    e        = orbits[:,1].copy()
    incl     = orbits[:,2].copy()
    longnode = orbits[:,3].copy()
    argperi  = orbits[:,4].copy()
    tperi    = orbits[:,5].copy()
    
    a = q/(1-e)
    meanmotion = np.sqrt(GMsun/(a*a*a), dtype=np.double)
    meananom = meanmotion*(jd_tt - tperi - lt)

    positions, velocities = kc.cartesian_vectors(orbits.shape[0], GMsun, a, e, incl, longnode, argperi, meananom)
    return np.array(positions).reshape((-1,3)), np.array(velocities).reshape((-1,3))

def rotate_matrix(ecl):
    ce = np.cos(ecl)
    se = np.sin(-ecl)
    rotmat = np.array([[1.0, 0.0, 0.0],
                  [0.0,  ce,  se],
                  [0.0, -se,  ce]])
    return rotmat

# This routine gets the geocentric position vectors of all
# the minor planets, corrected for light time.  It assumes
# that 3 iterations is adequate.  That assumption needs to 
# be checked.
rot_mat = rotate_matrix(Constants.ecl)
def getTopocentricPositions(GMsun, orbits, jd_utc, obsCode=None):
    jd_tt = EOP.jdTT(jd_utc)
    # Assume that the light time is 20 minutes as a first guess.    
    lt = np.ones(orbits.shape[0])*(20./(24.*60.))
    pos = Observatories.getObservatoryPosition(obsCode, jd_utc)
    lt_prev = np.zeros(orbits.shape[0])*(20./(24.*60.))
    i = 0
    while np.max(np.abs(lt-lt_prev))>1e-8:
        lt_prev = lt.copy()
        positions, _ = getPVs(GMsun, orbits, jd_tt, lt)
        heliopos = np.dot(rot_mat, positions.T).T
        geopos = heliopos-pos
        distances = np.linalg.norm(geopos, axis=1)
        lt = distances/Constants.speed_of_light
        i += 1
        #print i, lt-lt_prev
    return geopos, distances, heliopos

class PixelsOrbits():
    def __init__(self, pixel_path='pixel_dicts_test', orb_path='orb_dicts_test', mpc=True):

        self.orb_dict={}        
        self.pixel_dict = {}

        if mpc:
            mpc_files = sorted([join(orb_path, f) for f in listdir(orb_path) if isfile(join(orb_path, f)) & f.endswith('.txt')])
            #mpc_files = sorted([join(orb_path, f) for f in listdir(orb_path) if isfile(join(orb_path, f)) & f.endswith('.hdf5')])
            orb_epochs = sorted([f.rstrip('.txt').lstrip('mpn_') for f in listdir(orb_path) if isfile(join(orb_path, f)) & f.endswith('.txt')])
            orb_epochs_jd = np.array([MPC.yrmndy2JD(MPC.convertEpoch(epoch))  for epoch in orb_epochs])

            for i, (orb_epoch_jd, mpc_file) in enumerate(zip(orb_epochs_jd, mpc_files)):
                print  i, orb_epoch_jd, mpc_file
                self.orb_dict[orb_epoch_jd] = readMPCfile(filename=mpc_file)
                #self.orb_dict[orb_epoch_jd] = readMPC_HDF5_file(filename=mpc_file)

            print "here"
            pixel_files = [join(pixel_path, f) for f in listdir(pixel_path) if isfile(join(pixel_path, f)) & f.endswith('.5')]
            self.pixel_epochs_jd = sorted([float(pixel_file.split('_')[-1]) for pixel_file in pixel_files])

            ## Build a dictionary of the pixel dictionaries.
            for i, pixel_file in enumerate(pixel_files):
                pixel_epoch_jd = float(pixel_file.split('_')[-1])
                print i, pixel_epoch_jd, pixel_file
                self.pixel_dict[pixel_epoch_jd] = readHEALPix(pixel_file)

        else:
            mpc_files = sorted([join(orb_path, f) for f in listdir(orb_path) if isfile(join(orb_path, f)) & f.endswith('.pickle')])
            orb_epochs = sorted([f.rstrip('.txt').lstrip('mpn_') for f in listdir(orb_path) if isfile(join(orb_path, f)) & f.endswith('.pickle')])
            orb_epochs_jd = np.array([MPC.yrmndy2JD(MPC.convertEpoch(epoch))  for epoch in orb_epochs])

            for i, (orb_epoch_jd, mpc_file) in enumerate(zip(orb_epochs_jd, mpc_files)):
                print  i, orb_epoch_jd, mpc_file
                self.orb_dict[orb_epoch_jd] = readMPCorbits(filename=mpc_file)

            pixel_files = [join(pixel_path, f) for f in listdir(pixel_path) if isfile(join(pixel_path, f)) & f.endswith('.pickle')]
            self.pixel_epochs_jd = sorted([float(pixel_file.split('_')[-2]) for pixel_file in pixel_files])

            ## Build a dictionary of the pixel dictionaries.
            for i, pixel_file in enumerate(pixel_files):
                pixel_epoch_jd = float(pixel_file.split('_')[-2])
                print i, pixel_epoch_jd, pixel_file
                self.pixel_dict[pixel_epoch_jd] = readPixels(pixel_file)


    def getOrbitDict(self, jd_utc):
        closest_jd = min(self.orb_dict.keys(), key=lambda x:abs(x-jd_utc))
        return self.orb_dict[closest_jd]

    def getPixelDict(self, jd_utc):
        closest_jd = min(self.pixel_dict.keys(), key=lambda x:abs(x-jd_utc))
        return self.pixel_dict[closest_jd]
    

pixelsOrbits = PixelsOrbits()

def MPChecker(jd_utc, ra, dec, sr, pixelsOrbits=pixelsOrbits, nside=32, nested=True, obsCode=None):
    
    orb_dict   = pixelsOrbits.getOrbitDict(jd_utc)
    pixel_dict = pixelsOrbits.getPixelDict(jd_utc)
    
    # Find the healpix regions covered by the search
    # Probably should expand the search radius to accommodate 
    # daily motion.
    phi   = ra*math.pi/180.
    theta = math.pi/2.0 - dec*math.pi/180.
    if sr<1.0:
        search_radius = 1.0*math.pi/180.
    else:
        search_radius = sr*math.pi/180.
    vec    = hp.ang2vec(theta, phi)
    pixels = hp.query_disc(nside, vec, search_radius, nest=nested, inclusive=True)

    # Find the unique set of minor planets in the healpix regions
    mp_set = set()
    for pix in pixels:
        mps = pixel_dict[pix]
        mp_set.update(mps)

    orbits = np.zeros((len(mp_set), 6))
    orbit_table = {}
    H = np.zeros(len(mp_set))
    G = np.zeros(len(mp_set))
    for i, mp in enumerate(mp_set):
        orbit_table[mp] = i
        H[i], G[i], orbits[i] = orb_dict[mp]
        orbits[i] = np.array(orbits[i])

    # Compute the positions of those minor planets
    geopos, distances, heliopos = getTopocentricPositions(Constants.GMsun, orbits, jd_utc, obsCode=obsCode)
    
    # Now determine which MPs are actually in the search region
    geopos_normed=geopos/distances[:, np.newaxis]
    angles = np.degrees(np.arccos(np.dot(geopos_normed, vec)))

    r = np.sqrt(heliopos[:,0]*heliopos[:,0] + heliopos[:,1]*heliopos[:,1] + heliopos[:,2]*heliopos[:,2])
    d = np.sqrt(geopos[:,0]*geopos[:,0] + geopos[:,1]*geopos[:,1] + geopos[:,2]*geopos[:,2])

    phase_angle = 180./np.pi * np.arccos( (heliopos[:,0]*geopos[:,0] + heliopos[:,1]*geopos[:,1] + heliopos[:,2]*geopos[:,2])/(r*d))

    mag = MPC.H_alpha(H, G, phase_angle) + 5.0*np.log10(r*d)

    desig_result = []
    geopos_result = []
    distances_result = []
    angles_result = []
    mag_result = []
    alpha_result = []

    for mp in mp_set:
        i = orbit_table[mp]
        if angles[i]<sr:
            desig_result.append(mp)
            geopos_result.append(geopos[i])
            distances_result.append(distances[i])
            angles_result.append(angles[i])
            mag_result.append(mag[i])
            alpha_result.append(phase_angle[i])
    
    return desig_result, geopos_result, distances_result, angles_result, mag_result, alpha_result

if __name__ == '__main__':

  # Get arguments from command line and do checks
  try:
      inlist = sys.argv[1]
  except:
    print 'MPChecker inlist'
    sys.exit()

fov = 1.7 # search radius in degrees

with open(inlist, "r") as file:
    start = time.time()
    for line in file:
        pieces = line.split()
        if pieces[0] == 'Date':
            continue

        smf = pieces[2]
        jd_utc = float(pieces[6])
        ra = float(pieces[12])
        dec = float(pieces[13])

        desigs, geopos, distances, angles, mags, alphas = MPChecker(jd_utc, ra, dec, fov, obsCode='F51')
        for desig, pos, dist, angle, mag, alpha in zip(desigs, geopos, distances, angles, mags, alphas):
            x, y, z = pos
            x /= dist
            y /= dist
            z /= dist
            ra = np.degrees(np.arctan2(y, x))
            dec = np.degrees(np.arcsin(z))
            print desig, ra, dec, mag, angle
    end = time.time()
    print end - start






