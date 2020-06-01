# Retrives basic Gaia DR2 info given a source_id, coordinates, or Simbad name.

__version__ = '1.0'
__author__ = 'Ingrid Pelisoli'

# Importing relevant packages:

import numpy as np
import sys
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, Distance, Angle
from astropy.time import Time
from astroquery.mast import Observations
from astroquery.gaia import Gaia
from astropy.io.votable import parse_single_table

def plot_cmd(MG, bprp):
    table = parse_single_table("../TESS-LS/SampleC.vot")
    data = table.array

    s_MG = 5 + 5*np.log10(table.array['parallax']/1000) + table.array['phot_g_mean_mag']
    s_bprp = table.array['bp_rp']

    plt.scatter(s_bprp,s_MG,c='0.75', s=0.5, zorder=0)
    plt.gca().invert_yaxis()
    plt.title('$Gaia$ HR-diagram')
    plt.plot(bprp,MG,'or',markersize=10,zorder=2)
    plt.ylabel('$M_G$')
    plt.xlabel('$G_{BP}-G_{RP}$')
    plt.show()

def find_in_gaia(ra, dec):
    # First do a large search using 30 arcsec
    coord = SkyCoord(ra=ra, dec=dec,
                     unit=(u.degree, u.degree), frame='icrs')
    radius = u.Quantity(30.0, u.arcsec)
    q = Gaia.cone_search_async(coord, radius)
    gaia = q.get_results()
    gaia = gaia[ np.nan_to_num(gaia['parallax']) > 0 ]
    warning = (len(gaia) == 0)

    # Then propagate the Gaia coordinates to 2000, and find the best match to the
    # input coordinates
    if not warning:
        ra2015 = np.array(gaia['ra']) * u.deg
        dec2015 = np.array(gaia['dec']) * u.deg
        parallax = np.array(gaia['parallax']) * u.mas
        pmra = np.array(gaia['pmra']) * u.mas/u.yr
        pmdec = np.array(gaia['pmdec']) * u.mas/u.yr
        c2015 = SkyCoord(ra=ra2015, dec=dec2015,
                         distance=Distance(parallax=parallax, allow_negative=True),
                         pm_ra_cosdec=pmra, pm_dec=pmdec,
                         obstime=Time(2015.5, format='decimalyear'))
        c2000 = c2015.apply_space_motion(dt=-15.5 * u.year)

        idx, sep, _ = coord.match_to_catalog_sky(c2000)

        # The best match object
        best = gaia[idx]
        gaia_id = best['source_id']

        MG = 5 + 5*np.log10(best['parallax']/1000) + best['phot_g_mean_mag']
        bprp = best['bp_rp']

        gaia_id = np.int(gaia_id)
        G = np.float(best['phot_g_mean_mag'])
        MG = np.float(MG)
        bprp = np.float(bprp)

        return gaia_id, G, MG, bprp

if (len(sys.argv) == 3):
    ra = np.float(sys.argv[1])
    dec = np.float(sys.argv[2])
    gaia_id, G, MG, bprp = find_in_gaia(ra, dec)
    print("source_id            RA2015.5 DEC2015.5 G     MG    bp_rp")
    print("%20s %8.3f %9.3f %5.2f %5.2f %6.3f" % (gaia_id, ra, dec, G, MG, bprp))
    plot_cmd(MG, bprp)
elif (len(sys.argv) == 2):
    try:
        gaia_id = np.int(sys.argv[1])
        query = "select g.ra, g.dec, g.parallax, g.phot_g_mean_mag, g.bp_rp from gaiadr2.gaia_source as g where source_id = " + str(gaia_id)
        job = Gaia.launch_job(query=query)
        results = job.get_results()
        ra = np.float(results['ra'])
        dec = np.float(results['dec'])
        G = np.float(results['phot_g_mean_mag'])
        bprp = np.float(results['bp_rp'])
        MG = 5 + 5*np.log10(results['parallax']/1000) + results['phot_g_mean_mag']
        MG = np.float(MG)
        print("source_id            RA2015.5 DEC2015.5 G     MG    bp_rp")
        print("%20s %8.3f %9.3f %5.2f %5.2f %6.3f" % (gaia_id, ra, dec, G, MG, bprp))
        plot_cmd(MG, bprp)
    except ValueError:
        simbad_name = np.str(sys.argv[1])
        c = SkyCoord.from_name(simbad_name, frame='icrs')
        ra = c.ra.deg
        dec = c.dec.deg
        gaia_id, G, MG, bprp = find_in_gaia(ra, dec)
        print("source_id            RA2015.5 DEC2015.5 G     MG    bp_rp")
        print("%20s %8.3f %9.3f %5.2f %5.2f %6.3f" % (gaia_id, ra, dec, G, MG, bprp))
        plot_cmd(MG, bprp)
else:
    print("I don't understand your input.\n")
    print("Check README file for correct format.")
    sys.exit()
