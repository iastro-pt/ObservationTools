# -*- coding: utf-8 -*-
# from __future__ import print_function
import sys
import numpy as np
import datetime as dt
from dateutil import tz
from PyAstronomy import pyasl
from astropy.coordinates import SkyCoord
from astropy.coordinates import name_resolve
import ephem
import argparse

def _parser():
    parser = argparse.ArgumentParser(description='Plot altitudes of objects'
                                                 ' against time for a specific night')
    parser.add_argument('targets', help='E.g. HD20010 or HD20010,HD41248', nargs='+')
    parser.add_argument('-d', '--date', default='today',
                        help='Date in format YYYY-MM-DD (or YYYY if starobs). Default is today.')
    parser.add_argument('-s', '--site', default='esolasilla',
                        help='Observatory. Default is ESO La Silla. '
                             'Common codes are esoparanal, lapalma, keck, lco, Palomar, etc')
    parser.add_argument('-c', default=False, action='store_true',
                        help='Just print "target RA DEC" (to use in STARALT)')
    parser.add_argument('-m', '--mode', choices=['staralt', 'starobs'], default='staralt',
                        help='staralt: plot altitude against time for a particular night; starobs: plot how altitude changes over a year')
    return parser.parse_args()


def decdeg2dms(dd):
    """ Convert decimal degrees to deg,min,sec """
    is_positive = dd >= 0
    dd = abs(dd)
    minutes,seconds = divmod(dd*3600,60)
    degrees,minutes = divmod(minutes,60)
    degrees = degrees if is_positive else -degrees
    return (degrees,minutes,seconds)


def StarObsPlot(year=None, targets=None, observatory=None, print2file=False):
  """
    Plot the visibility of target.

    Parameters
    ----------
    year: int
        The year for which to calculate the visibility.
    targets: list
        List of targets.
        Each target should be a dictionary with keys 'name' and 'coord'.
        The key 'name' is a string, 'coord' is a SkyCoord object.
    observatory: string
        Name of the observatory that pyasl.observatory can resolve.
        Basically, any of pyasl.listObservatories().keys()
    print2file : boolean, optional
        If True, the plot will be dumped to a png-file.
  """

  import matplotlib.pylab as plt
  from mpl_toolkits.axes_grid1 import host_subplot
  from matplotlib.ticker import MultipleLocator
  from matplotlib.font_manager import FontProperties
  from matplotlib import rcParams
  rcParams['xtick.major.pad'] = 12
  font0 = FontProperties()
  font1 = font0.copy()
  font0.set_family('sans-serif')
  font0.set_weight('light')
  font1.set_family('sans-serif')
  font1.set_weight('medium')

  # set the observatory
  obs = pyasl.observatory(observatory)


  fig = plt.figure(figsize=(15,10))
  fig.subplots_adjust(left=0.07, right=0.8, bottom=0.15, top=0.88)
  ax = host_subplot(111)

  for n, target in enumerate(targets):

    target_coord = target['coord']
    target_ra = target_coord.ra.deg
    target_dec = target_coord.dec.deg

    jdbinsize = 1 # every day
    jd_start = pyasl.jdcnv(dt.datetime(year, 1, 1))
    jd_end = pyasl.jdcnv(dt.datetime(year, 12, 31))
    each_day = np.arange(jd_start, jd_end, jdbinsize)
    jds = []

    ## calculate the mid-dark times
    sun = ephem.Sun()
    for day in each_day:
      date_formatted = '/'.join([str(i) for i in pyasl.daycnv(day)[:-1]])
      s = ephem.Observer();
      s.date = date_formatted;
      s.lat = ':'.join([str(i) for i in decdeg2dms(obs['latitude'])])
      s.lon = ':'.join([str(i) for i in decdeg2dms(obs['longitude'])])
      jds.append(ephem.julian_date(s.next_antitransit(sun)))
    jds = np.array(jds)

    # Get JD floating point
    jdsub = jds - np.floor(jds[0])

    # Get alt/az of object
    altaz = pyasl.eq2hor(jds, np.ones_like(jds)*target_ra, np.ones_like(jds)*target_dec, \
                        lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])
    ax.plot( jdsub, altaz[0], '-', color='k')


    # label for each target
    plabel = "[{0:2d}]  {1!s}".format(n+1, target['name'])

    # number of target at the top of the curve
    altmax = np.argmax(altaz[0])
    ax.text(jdsub[altmax], altaz[0][altmax], str(n+1), color="b", fontsize=14, \
             fontproperties=font1, va="bottom", ha="center")

    if n+1 == 29:
      # too many?
      ax.text(1.1, 1.0-float(n+1)*0.04, "too many targets", ha="left", va="top", transform=ax.transAxes, \
              fontsize=10, fontproperties=font0, color="r")
    else:
      ax.text(1.1, 1.0-float(n+1)*0.04, plabel, ha="left", va="top", transform=ax.transAxes, \
              fontsize=12, fontproperties=font0, color="b")


  ax.text(1.1, 1.03, "List of targets", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")


  axrange = ax.get_xlim()
  months = range(1, 13)
  import calendar
  ndays = [0] + [calendar.monthrange(date, m)[1] for m in months]

  if axrange[1]-axrange[0] <= 1.0:
    jdhours = np.arange(0,3,1.0/24.)
    utchours = (np.arange(0,72,dtype=int)+12)%24
  else:
    jdhours = np.arange(0,3,1.0/12.)
    utchours = (np.arange(0,72, 2, dtype=int)+12)%24

  ax.set_xlim([0, 366])
  ax.set_xticks(np.cumsum(ndays)[:-1] + (np.array(ndays)/2.)[1:])
  ax.set_xticklabels(map(calendar.month_abbr.__getitem__, months), fontsize=10)

  # Make ax2 responsible for "top" axis and "right" axis
  ax2 = ax.twin()
  # Set upper x ticks
  ax2.set_xticks(np.cumsum(ndays))
  ax2.set_xlabel("Day")

  # Horizon angle for airmass
  airmass_ang = np.arange(5.,90.,5.)
  geo_airmass = pyasl.airmass.airmassPP(90.-airmass_ang)
  ax2.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(geo_airmass.size):
    airmassformat.append("{0:2.2f}".format(geo_airmass[t]))
  ax2.set_yticklabels(airmassformat, rotation=90)
  ax2.set_ylabel("Relative airmass", labelpad=32)
  ax2.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.015,-0.04, "Plane-parallel", transform=ax.transAxes, ha='left', \
           va='top', fontsize=10, rotation=90)

  ax22 = ax.twin()
  ax22.set_xticklabels([])
  ax22.set_frame_on(True)
  ax22.patch.set_visible(False)
  ax22.yaxis.set_ticks_position('right')
  ax22.yaxis.set_label_position('right')
  ax22.spines['right'].set_position(('outward', 25))
  ax22.spines['right'].set_color('k')
  ax22.spines['right'].set_visible(True)
  airmass2 = list(map(lambda ang: pyasl.airmass.airmassSpherical(90. - ang, obs['altitude']), airmass_ang))
  ax22.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(len(airmass2)): airmassformat.append("{0:2.2f}".format(airmass2[t]))
  ax22.set_yticklabels(airmassformat, rotation=90)
  ax22.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.045,-0.04, "Spherical+Alt", transform=ax.transAxes, ha='left', va='top', \
           fontsize=10, rotation=90)


  ax.set_ylim([0, 91])
  ax.yaxis.set_major_locator(MultipleLocator(15))
  ax.yaxis.set_minor_locator(MultipleLocator(5))
  yticks = ax.get_yticks()
  ytickformat = []
  for t in range(yticks.size):
    ytickformat.append(str(int(yticks[t]))+r"$^\circ$")
  ax.set_yticklabels(ytickformat, fontsize=16)
  ax.set_ylabel("Altitude", fontsize=18)
  yticksminor = ax.get_yticks(minor=True)
  ymind = np.where( yticksminor % 15. != 0. )[0]
  yticksminor = yticksminor[ymind]
  ax.set_yticks(yticksminor, minor=True)
  m_ytickformat = []
  for t in range(yticksminor.size):
    m_ytickformat.append(str(int(yticksminor[t]))+r"$^\circ$")
  ax.set_yticklabels(m_ytickformat, minor=True)
  ax.set_ylim([0, 91])

  ax.yaxis.grid(color='gray', linestyle='dashed')
  ax.yaxis.grid(color='gray', which="minor", linestyle='dotted')
  ax2.xaxis.grid(color='gray', linestyle='dotted')

  plt.text(0.5,0.95,"Visibility over {0!s}\n - altitudes at mid-dark time -".format(date), \
           transform=fig.transFigure, ha='center', va='bottom', fontsize=16)


  obsco = "Obs coord.: {0:8.4f}$^\circ$, {1:8.4f}$^\circ$, {2:4f} m".format(obs['longitude'], obs['latitude'], obs['altitude'])

  plt.text(0.01,0.97, obsco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  plt.text(0.01,0.95, obs['name'], transform=fig.transFigure, ha='left', va='center', fontsize=10)

  if print2file:
    pass
    # plt.savefig(outfile, format="png", dpi=300)
  else:
    plt.show()



def VisibilityPlot(date=None, targets=None, observatory=None, plotLegend=True, showMoonDist=True, print2file=False):
  """
    Plot the visibility of target.

    Parameters
    ----------
    date: datetime
        The date for which to calculate the visibility.
    targets: list
        List of targets.
        Each target should be a dictionary with keys 'name' and 'coord'.
        The key 'name' is aa string, 'coord' is a SkyCoord object.
    observatory: string
        Name of the observatory that pyasl.observatory can resolve.
        Basically, any of pyasl.listObservatories().keys()
    plotLegend: boolean, optional
        If True (default), show a legend.
    showMoonDist : boolean, optional
        If True (default), the Moon distance will be shown.
    print2file : boolean, optional
        If True, the plot will be dumped to a png-file.
  """

  try:
    import matplotlib
    import matplotlib.pylab as plt
    from mpl_toolkits.axes_grid1 import host_subplot
    from matplotlib.ticker import MultipleLocator
    from matplotlib.font_manager import FontProperties
    from matplotlib import rcParams
  except ImportError:
    print('matplotlib is not installed?')
    sys.exit(1)

  rcParams['xtick.major.pad'] = 12


  obs = pyasl.observatory(observatory)

  # observer = ephem.Observer()
  # observer.pressure = 0
  # observer.horizon = '-0:34'
  # observer.lat, observer.lon = obs['latitude'], obs['longitude']
  # observer.date = date
  # print(observer.date)
  # print(observer.previous_rising(ephem.Sun()))
  # print(observer.next_setting(ephem.Sun()))
  # print(observer.previous_rising(ephem.Moon()))
  # print(observer.next_setting(ephem.Moon()))
  # observer.horizon = '-6'
  # noon = observer.next_transit(ephem.Sun())
  # print(noon)
  # print(observer.previous_rising(ephem.Sun(), start=noon, use_center=True))
  # print()


  fig = plt.figure(figsize=(15,10))
  fig.subplots_adjust(left=0.07, right=0.8, bottom=0.15, top=0.88)
  ax = host_subplot(111)

  font0 = FontProperties()
  font1 = font0.copy()
  font0.set_family('sans-serif')
  font0.set_weight('light')
  font1.set_family('sans-serif')
  font1.set_weight('medium')


  for n, target in enumerate(targets):

    target_coord = target['coord']
    target_ra = target_coord.ra.deg
    target_dec = target_coord.dec.deg

    # JD array
    jdbinsize = 1.0/24./20.
    # jds = np.arange(allData[n]["Obs jd"][0], allData[n]["Obs jd"][2], jdbinsize)
    jd = pyasl.jdcnv(date)
    jd_start = pyasl.jdcnv(date)-0.5
    jd_end = pyasl.jdcnv(date)+0.5
    jds = np.arange(jd_start, jd_end, jdbinsize)
    # Get JD floating point
    jdsub = jds - np.floor(jds[0])
    # Get alt/az of object
    altaz = pyasl.eq2hor(jds, np.ones(jds.size)*target_ra, np.ones(jds.size)*target_dec, \
                        lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])
    # Get alt/az of Sun
    sun_position = pyasl.sunpos(jd)
    sun_ra, sun_dec = sun_position[1], sun_position[2]
    sunpos_altaz = pyasl.eq2hor(jds, np.ones(jds.size)*sun_ra, np.ones(jds.size)*sun_dec, \
                                lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])

    # Define plot label
    plabel = "[{0:2d}]  {1!s}".format(n+1, target['name'])

    # Find periods of: day, twilight, and night
    day = np.where( sunpos_altaz[0] >= 0. )[0]
    twi = np.where( np.logical_and(sunpos_altaz[0] > -18., sunpos_altaz[0] < 0.) )[0]
    night = np.where( sunpos_altaz[0] <= -18. )[0]

    if (len(day) == 0) and (len(twi) == 0) and (len(night) == 0):
      print
      print("VisibilityPlot - no points to draw")
      print

    mpos = pyasl.moonpos(jds)
    # mpha = pyasl.moonphase(jds)
    # mpos_altaz = pyasl.eq2hor(jds, mpos[0], mpos[1],
    #                            lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])
    # moonind = np.where( mpos_altaz[0] > 0. )[0]

    if showMoonDist:
      mdist = pyasl.getAngDist(mpos[0], mpos[1], np.ones(jds.size)*target_ra, \
                              np.ones(jds.size)*target_dec)
      bindist = int((2.0/24.)/jdbinsize)
      firstbin = np.random.randint(0,bindist)
      for mp in range(0, int(len(jds)/bindist)):
        bind = firstbin+mp*bindist
        if altaz[0][bind]-1. < 5.: continue
        ax.text(jdsub[bind], altaz[0][bind]-1., str(int(mdist[bind]))+r"$^\circ$", ha="center", va="top", \
                fontsize=8, stretch='ultra-condensed', fontproperties=font0, alpha=1.)


    if len(twi) > 1:
      # There are points in twilight
      linebreak = np.where( (jdsub[twi][1:]-jdsub[twi][:-1]) > 2.0*jdbinsize)[0]
      if len(linebreak) > 0:
        plotrjd = np.insert(jdsub[twi], linebreak+1, np.nan)
        plotdat = np.insert(altaz[0][twi], linebreak+1, np.nan)
        ax.plot( plotrjd, plotdat, "-", color='#BEBEBE', linewidth=1.5)
      else:
        ax.plot( jdsub[twi], altaz[0][twi], "-", color='#BEBEBE', linewidth=1.5)

    ax.plot( jdsub[night], altaz[0][night], '.k', label=plabel)
    ax.plot( jdsub[day], altaz[0][day], '.', color='#FDB813')

    altmax = np.argmax(altaz[0])
    ax.text( jdsub[altmax], altaz[0][altmax], str(n+1), color="b", fontsize=14, \
             fontproperties=font1, va="bottom", ha="center")

    if n+1 == 29:
      ax.text( 1.1, 1.0-float(n+1)*0.04, "too many targets", ha="left", va="top", transform=ax.transAxes, \
              fontsize=10, fontproperties=font0, color="r")
    else:
      ax.text( 1.1, 1.0-float(n+1)*0.04, plabel, ha="left", va="top", transform=ax.transAxes, \
              fontsize=12, fontproperties=font0, color="b")

  ax.text( 1.1, 1.03, "List of targets", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")

  axrange = ax.get_xlim()
  ax.set_xlabel("UT [hours]")

  if axrange[1]-axrange[0] <= 1.0:
    jdhours = np.arange(0,3,1.0/24.)
    utchours = (np.arange(0,72,dtype=int)+12)%24
  else:
    jdhours = np.arange(0,3,1.0/12.)
    utchours = (np.arange(0,72, 2, dtype=int)+12)%24
  ax.set_xticks(jdhours)
  ax.set_xlim(axrange)
  ax.set_xticklabels(utchours, fontsize=18)

  # Make ax2 responsible for "top" axis and "right" axis
  ax2 = ax.twin()
  # Set upper x ticks
  ax2.set_xticks(jdhours)
  ax2.set_xticklabels(utchours, fontsize=18)
  ax2.set_xlabel("UT [hours]")

  # Horizon angle for airmass
  airmass_ang = np.arange(5.,90.,5.)
  geo_airmass = pyasl.airmass.airmassPP(90.-airmass_ang)
  ax2.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(geo_airmass.size):
    airmassformat.append("{0:2.2f}".format(geo_airmass[t]))
  ax2.set_yticklabels(airmassformat, rotation=90)
  ax2.set_ylabel("Relative airmass", labelpad=32)
  ax2.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.015,-0.04, "Plane-parallel", transform=ax.transAxes, ha='left', \
           va='top', fontsize=10, rotation=90)

  ax22 = ax.twin()
  ax22.set_xticklabels([])
  ax22.set_frame_on(True)
  ax22.patch.set_visible(False)
  ax22.yaxis.set_ticks_position('right')
  ax22.yaxis.set_label_position('right')
  ax22.spines['right'].set_position(('outward', 25))
  ax22.spines['right'].set_color('k')
  ax22.spines['right'].set_visible(True)
  airmass2 = list(map(lambda ang: pyasl.airmass.airmassSpherical(90. - ang, obs['altitude']), airmass_ang))
  ax22.set_yticks(airmass_ang)
  airmassformat = []
  for t in airmass2:
      airmassformat.append("{0:2.2f}".format(t))
  ax22.set_yticklabels(airmassformat, rotation=90)
  ax22.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.045,-0.04, "Spherical+Alt", transform=ax.transAxes, ha='left', va='top', \
           fontsize=10, rotation=90)

  ax3 = ax.twiny()
  ax3.set_frame_on(True)
  ax3.patch.set_visible(False)
  ax3.xaxis.set_ticks_position('bottom')
  ax3.xaxis.set_label_position('bottom')
  ax3.spines['bottom'].set_position(('outward', 50))
  ax3.spines['bottom'].set_color('k')
  ax3.spines['bottom'].set_visible(True)

  ltime, ldiff = pyasl.localtime.localTime(utchours, np.repeat(obs['longitude'], len(utchours)))
  jdltime = jdhours - ldiff/24.
  ax3.set_xticks(jdltime)
  ax3.set_xticklabels(utchours)
  ax3.set_xlim([axrange[0],axrange[1]])
  ax3.set_xlabel("Local time [hours]")

  ax.set_ylim([0, 91])
  ax.yaxis.set_major_locator(MultipleLocator(15))
  ax.yaxis.set_minor_locator(MultipleLocator(5))
  yticks = ax.get_yticks()
  ytickformat = []
  for t in range(yticks.size): ytickformat.append(str(int(yticks[t]))+r"$^\circ$")
  ax.set_yticklabels(ytickformat, fontsize=16)
  ax.set_ylabel("Altitude", fontsize=18)
  yticksminor = ax.get_yticks(minor=True)
  ymind = np.where( yticksminor % 15. != 0. )[0]
  yticksminor = yticksminor[ymind]
  ax.set_yticks(yticksminor, minor=True)
  m_ytickformat = []
  for t in range(yticksminor.size): m_ytickformat.append(str(int(yticksminor[t]))+r"$^\circ$")
  ax.set_yticklabels(m_ytickformat, minor=True)
  ax.set_ylim([0, 91])

  ax.yaxis.grid(color='gray', linestyle='dashed')
  ax.yaxis.grid(color='gray', which="minor", linestyle='dotted')
  ax2.xaxis.grid(color='gray', linestyle='dotted')

  plt.text(0.5,0.95,"Visibility on {0!s}".format(date.date()), \
           transform=fig.transFigure, ha='center', va='bottom', fontsize=20)

  if plotLegend:
    line1 = matplotlib.lines.Line2D((0,0),(1,1), color='#FDB813', linestyle="-", linewidth=2)
    line2 = matplotlib.lines.Line2D((0,0),(1,1), color='#BEBEBE', linestyle="-", linewidth=2)
    line3 = matplotlib.lines.Line2D((0,0),(1,1), color='k', linestyle="-", linewidth=2)

    lgd2 = plt.legend((line1,line2,line3),("day","twilight","night",), \
                        bbox_to_anchor=(0.88, 0.13), loc='best', borderaxespad=0.,prop={'size':12}, fancybox=True)
    lgd2.get_frame().set_alpha(.5)

  obsco = "Obs coord.: {0:8.4f}$^\circ$, {1:8.4f}$^\circ$, {2:4f} m".format(obs['longitude'], obs['latitude'], obs['altitude'])

  plt.text(0.01,0.97, obsco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  plt.text(0.01,0.95, obs['name'], transform=fig.transFigure, ha='left', va='center', fontsize=10)

  if print2file:
    pass
    # plt.savefig(outfile, format="png", dpi=300)
  else:
    plt.show()



if __name__ == '__main__':
  args = _parser()

  target_names = args.targets[0].split(',')
  # print(target_names)

  ## Get coordinates for all the targets
  targets = []
  for target_name in target_names:
    try:
      targets.append({'name': target_name, 'coord': SkyCoord.from_name(target_name)})
    except name_resolve.NameResolveError as e:
      print('Could not find target: {0!s}'.format(target_name))

  ## Just print coordinates in STARALT format and exit
  if args.c:
    print('Coordinates for {0!s}\n'.format(args.targets[0]))
    for target in targets:
      ## name hh mm ss ±dd mm ss
      out = '{0!s}'.format(target['name'])
      ra = target['coord'].ra.hms
      out += ' {0:02d} {1:02d} {2:5.3f}'.format(ra.h, ra.m, ra.s)
      dec = target['coord'].dec.dms
      out += ' {0:02d} {1:02d} {2:5.3f}'.format(dec.d, dec.m, dec.s)
      print(out)

    sys.exit(0)

  ## Actually calculate the visibility curves
  print('Calculating visibility for {0!s}'.format(args.targets[0]))

  if args.date == 'today':
    if args.mode == 'staralt':
      today = dt.datetime.now() # now() gives the current time which we don't want
      date = dt.datetime(today.year, today.month, today.day, tzinfo=tz.tzutc())
      print(date.date())
    elif args.mode == 'starobs':
      date = dt.datetime.now().year
      print(date)
  else:
    if args.mode == 'staralt':
      if "-" not in args.date:
        raise ValueError("Date needs to be provided as YYYY-MM-DD for staralt mode.")
      ymd = [int(i) for i in args.date.split('-')]
      date = dt.datetime(*ymd)
      print(date.date())
    elif args.mode == 'starobs':
      if "-" in args.date:
        date = int(args.date.split('-')[0])
      else:
        date = int(args.date)
      print(date)

  ## Find observatory
  available_sites = pyasl.listObservatories(show=False)
  if args.site not in available_sites.keys():
    print('"{0!s}" is not a valid observatory code. Try one of the following:\n'.format(args.site))

    maxCodeLen = max(map(len, available_sites.keys()))
    print(("{0:"+str(maxCodeLen)+"s}     ").format("Code") + "Observatory name")
    print("-" * (21+maxCodeLen))
    for k in sorted(available_sites.keys(), key=lambda s: s.lower()):
      print(("{0:"+str(maxCodeLen)+"s} --- ").format(k) + available_sites[k]["name"])

    sys.exit(1)


  print(pyasl.observatory(args.site)['name'])

  if args.mode == 'staralt':
    ## Plot visibility
    VisibilityPlot(date=date, targets=targets, observatory=args.site)
  elif args.mode == 'starobs':
    StarObsPlot(year=date, targets=targets, observatory=args.site)
