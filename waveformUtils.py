import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import math
import matplotlib as mpl
from matplotlib import mlab
from matplotlib.colors import Normalize
from copy import copy
import matplotlib.gridspec as gridspec
from obspy import read
from obspy import Stream
from obspy.clients.fdsn import Client
from pathlib import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plotRelativeAlign(s0, masterIndex=0, normalize=True, figsize=(8, 10), shifts=None, picks=None, polarity=False, handle=None):
    """
    Plot obspy stream object aligned based on a single waveform.

    Parameters
    ----------
    s0: Obspy Stream
        Stream containing events
    masterIndex : integer
        Integer of the waveform within the stream object to use as the master for cross-correlation.  Defaults to zero (the first trace).
    normalize : boolean
        Choose whether or not to normalize the stream
    figsize: tuple
        Figure size (width, length).  Default to (8,10).
    shifts: list
        List of shifts in seconds to apply to the data.  Must be same length as s0.
    picks: list of UTCDateTime
        List of UTCDateTime objects on same length as s0.  If no pick is available, set to None
    polarity: boolean
        Boolean on whether or not to plot polarity (red is positive, blue is negative, black is none).
        Must have "polarity" field in stats. (For example: Snow[0].stats.polarity = nowpick.polarity)
    handle: axis handle
        pyplot axis handle to put plot within

    Outputs
    ---------
    f: matplotlib figure handle
        Figure handle from matplotlib to adjust the plot if necessary
    shifts: list
        Shifts (in seconds) used to get the waveforms into alignment
    """
    import numpy as np
    import obspy.signal.cross_correlation as xc

    plotPicks = False
    if picks:
        plotPicks = True
        if len(picks) != len(s0):
            print ('Not enough picks for the number of traces.')
            return f, shifts
    s1 = s0.copy()
    if not handle:
        f = plt.figure(figsize=figsize)
        ax1 = f.add_subplot()
    else:
        ax1 = handle
        plt.sca(ax1)
    cntr = 0
    maxtime = 0
    if normalize:
        s1.normalize()
    if shifts:
        if len(shifts) != len(s0):
            print ('Not enough shifts for number of traces within stream.')
            return f, shifts
        else:
            if plotPicks:
                for st, shift, pick in zip(s1, shifts, picks):
                    t = st.times()
                    t = t - shift
                    if polarity and pick:
                        if st.stats.polarity:
                            if st.stats.polarity == 'positive':
                                plt.plot(t, st.data + cntr, '-r')
                            elif st.stats.polarity == 'negative':
                                plt.plot(t, st.data + cntr, '-b')
                        else:
                            plt.plot(t, st.data + cntr, '-k')
                    else:
                        plt.plot(t, st.data + cntr, '-k')
                    if pick:
                        relptime = pick-st.stats.starttime-shift
                        plt.plot(relptime, cntr, 'or')
                    cntr -= 1
            else:
                for st, shift in zip(s1, shifts):
                    t = st.times()
                    t = t - shift
                    try:
                        if polarity:
                            if st.stats.polarity:
                                if st.stats.polarity == 'positive':
                                    plt.plot(t, st.data + cntr, '-r')
                                elif st.stats.polarity == 'negative':
                                    plt.plot(t, st.data + cntr, '-b')
                            else:
                                plt.plot(t, st.data + cntr, '-k')
                        else:
                            plt.plot(t, st.data + cntr, '-k')
                    except:
                        plt.plot(t, st.data + cntr, '-k')
                    cntr -= 1
    else:
        shifts = []
        if plotPicks:
            for st,pick in zip(s1, picks):  # Loop over waveforms to get time shifts (if not provided)
                t = st.times()
                if max(t) > maxtime:
                    maxtime = max(t)
                cc = xc.correlate(st.data, s1[masterIndex].data, int(len(s1[masterIndex].data) / 2))
                shift, value = xc.xcorr_max(cc)
                t = t - (shift * st.stats.delta)
                shifts.append(shift * st.stats.delta)
                if polarity and pick:
                    if st.stats.polarity:
                        if st.stats.polarity == 'positive':
                            plt.plot(t, st.data + cntr, '-r')
                        elif st.stats.polarity == 'negative':
                            plt.plot(t, st.data + cntr, '-b')
                    else:
                        plt.plot(t, st.data + cntr, '-k')
                else:
                    plt.plot(t, st.data + cntr, '-k')
                if pick:
                    relptime = pick-st.stats.starttime-(shift*st.stats.delta)
                    plt.plot(relptime, cntr, 'or')
                cntr -= 1
        else:
            for st in s1:  # Loop over waveforms to get time shifts (if not provided)
                t = st.times()
                if max(t) > maxtime:
                    maxtime = max(t)
                cc = xc.correlate(st.data, s1[masterIndex].data, int(len(s1[masterIndex].data) / 2))
                shift, value = xc.xcorr_max(cc)
                t = t - (shift * st.stats.delta)
                shifts.append(shift * st.stats.delta)
                try:
                    if polarity:
                        if st.stats.polarity:
                            if st.stats.polarity == 'positive':
                                plt.plot(t, st.data + cntr, '-r')
                            elif st.stats.polarity == 'negative':
                                plt.plot(t, st.data + cntr, '-b')
                        else:
                            plt.plot(t, st.data + cntr, '-k')
                    else:
                        plt.plot(t, st.data + cntr, '-k')
                except:
                    plt.plot(t, st.data + cntr, '-k')
                cntr -= 1

    ts = [nowt.stats.starttime for nowt in s1]
    plt.yticks(np.arange(-len(s1) + 1, 1), [nowt.strftime('%Y%m%d_%H:%M:%S') for nowt in reversed(ts)])
    return plt.gcf(), shifts


def plotRelative(s0, normalize=True, figsize=(8, 10), handle=None, fill_between=False, reduceAmp=None, linewidth=1, ylabels='date'):
    """
    Plot obspy stream object aligned based on a single waveform.

    Parameters
    ----------
    s0: Obspy Stream
        Stream containing events
    normalize : boolean
        Choose whether or not to normalize the stream
    figsize: tuple
        Figure size (width, length).  Default to (8,10).
    handle: axis handle
        pyplot axis handle to put plot within
    fill_between: boolean
        determines whether or not to fill the area between the local zero and trace w/
        gray.  Good for envelopes.
    reduceAmp: int (1-100)
        factor with which to reduce the amplitude by.  default uses all of the space between the separate events, 
        which can lead to some overlap.  If, using, start at 75 and go down from there.
    Outputs
    ---------
    f: matplotlib figure handle
        Figure handle from matplotlib to adjust the plot if necessary
    """
    import numpy as np

    s1 = s0.copy()
    if not handle:
        f = plt.figure(figsize=figsize)
        ax1 = f.add_subplot()
    else:
        ax1 = handle
        plt.sca(ax1)
    cntr = 0
    maxtime = 0
    if normalize:
        s1.normalize()
    if reduceAmp:
        for tr in s1:
            tr.data = tr.data*(reduceAmp/100.0)
    for st in s1:  # Loop over waveforms to get time shifts
        t = st.times()
        if max(t) > maxtime:
            maxtime = max(t)
        ax1.plot(t, st.data + cntr, '-k', lw=linewidth)
        if fill_between:
            ax1.fill_between(t, st.data+cntr, y2=cntr, color='gray', alpha=0.7 )
        cntr -= 1

    if ylabels=='date':
        ts = [nowt.stats.starttime for nowt in s1]
        plt.yticks(np.arange(-len(s1) + 1, 1), [nowt.strftime('%Y%m%d_%H:%M:%S') for nowt in reversed(ts)])
    else:
        plt.yticks(np.arange(-len(s1) + 1, 1), [tr.id for tr in reversed(s0)])
    
    ax1.set_xlim([0, maxtime])
    ax1.set_xlabel('Relative Time, s')
    plt.tight_layout()
    return plt.gcf()


def interferogramPlot(stream, normalize=True, figsize=(8, 10), shift=False, masterIndex=0):
    """
    
    :param stream: 
    :param normalize: 
    :param figsize: 
    :param shift: 
    :param masterIndex: 
    :return: 
    """
    # Use imshow w/ data matrix (need to pad edges)
    # plt.imshow(data.T, extent=xextent, aspect="auto",
    #                    interpolation='none', origin='lower', cmap='seismic',
    #                    vmin=-vmax, vmax=vmax)


def seisfft(traceData):
    """
    Takes an obspy trace object and calculates a spectra.  Make sure that the
    trace is demeaned and tapered.
    
    Parameters
    ---------------
    traceData: Obspy trace object
        Obspy trace object, typically seismic data
        
    Outputs
    ---------------
    f: array
        array of frequencies of spectra
    g: array
        discrete fourier spectrum values (real and imaginary)
    
    
    To plot:
        import matplotlib.pyplot as plt
        plt.plot(f, abs(g))
    """
    from numpy import fft
    data = traceData.data
    stats = traceData.stats
    npts = stats['npts']
    delta = stats['delta']
    g = fft.rfft(data)
    f = fft.rfftfreq(npts, delta)
    return f, g

def stackSpectraContinuous(S, window_length=60, step=30):

    import waveformUtils
    initial = True
    for tr in S:
        for w_st in tr.slide(window_length=window_length, step=step):
            F1, G1 = waveformUtils.seisfft(w_st)
            if initial:
                ffts = np.array(abs(G1))
                initial = False
            else:
                ffts = np.vstack((ffts, abs(G1)))

    fftStack = np.sum(ffts, axis=0) / len(ffts)
    return F1, fftStack

def plotSpectraTime(S1, minfreq=1, maxfreq=10, winlength=600, step=600, 
                    spacing=2, ampScalar=3, normalize=True, specWin=60,
                    specStep=30, rsamWin=1, saveFig=False):
    '''
    # Plot stacked spectra with respect to time
    :param S1: Stream (single channel (use st.select()))
    :param minfreq: minimum frequency to plot (default=1)
    :param maxfreq: maximum frequency to plot (default=10)
    :param winlength: number of seconds to stack spectra over (default=600)
    :param step: number of seconds to step (default=300)
    :param specWin: number of seconds to use for spectra (default=60)
    :param specStep: number of seconds to step between each spectra (default=30)
    :param spacing: scalar to increase vertical spacing between traces (default=2)
    :param ampScalar: scalar to increase amplitude of individual traces (default=3)
    :param normalize: boolean. True means normalize each trace (between minFreq and maxFreq). False means normalize to a global maximum. (default=True)
    :param rsamWin: window length in minutes (default=1)
    :param saveFig: save eps figure in current directory (default=False)
    :return: f, spectraF, spectraG, sttimes, tvec, RSAM 
    '''
    
    spectraG = []
    spectraF = []
    sttimes = []
    # Calculate spectra in sliding windows
    for wst in S1.slide(window_length=winlength, step=step, include_partial_windows=True):
        F, G = stackSpectraContinuous(wst.detrend(type='demean'), window_length=specWin, step=specStep)
        spectraF.append(F[F >= minfreq])
        spectraG.append(np.abs(G[F >= minfreq]))
        sttimes.append(wst[0].stats.starttime)

    f = plt.figure(figsize=(15, 20))
    ax = plt.subplot2grid((1, 6), (0, 0), rowspan=1, colspan=4)
    ax2 = plt.subplot2grid((1, 6), (0, 4), rowspan=1, colspan=2)
    cntr = 0
    # Find global max of spectra
    globalmax = 0
    if not normalize:
        for gnow in spectraG:
            if max(gnow) > globalmax:
                globalmax = max(gnow)
    # Plot spectra
    for fnow, gnow in zip(spectraF, spectraG):
        # Create normalized y
        if normalize:
            Y = ampScalar * (gnow / max(gnow)) + cntr
        else:
            Y = ampScalar * (gnow / globalmax) + cntr
        ax.plot(fnow, Y, '-k')
        cntr -= spacing
    # plt.yticks(np.arange((cntr+2) + 1, 2), [nowt.strftime('%Y%m%d_%H:%M:%S') for nowt in reversed(sttimes)])
    ax.set_ylim(cntr+spacing, ampScalar)
    yt = ax.get_yticks()
    ylabels = []
    for nt in yt:
        nowt = S1[0].stats.starttime + (-1 * (nt / spacing) * step)
        ylabels.append(nowt.strftime('%Y%m%d_%H:%M:%S'))
    ax.set_yticklabels(ylabels)
    ax.set_xlabel('Frequency, HZ')
    ax.set_ylabel('Normalized Spectra Start Time, UTC')
    ax.set_xlim(minfreq, maxfreq)
    ax.grid()

    # Get RSAM
    S2 = S1.copy().detrend(type='demean').taper(max_percentage=None, max_length=5).filter(type='bandpass', freqmin=minfreq,
                                                                                          freqmax=maxfreq)
    tvec, rsam = RSAM(S2[0], 1, timePos='start')
    plotvec = [nowt.matplotlib_date for nowt in tvec]
    ax2.plot(rsam, plotvec, '-k')
    ax2.set_ylim(min(plotvec), max(plotvec))
    ax2.set_xlim(0, max(rsam))
    ax2.set_xlim(0, np.mean(rsam)*3)
    ax2.invert_yaxis()
    ax2.set_yticklabels([])
    ax2.set_xlabel('RSAM Amplitude')
    f.show()

    if saveFig:
        filename = 'spectraFig_{}_{}_{}.eps'.format(S2[0].stats.station, S2[0].stats.starttime.strftime('%Y%m%d%H%M%S'),
                                                S2[0].stats.endtime.strftime('%Y%m%d%H%M%S'))
        f.savefig(filename)
        
    return f, spectraF, spectraG, sttimes, tvec, RSAM
    
def RSAM(trace, window, timePos='middle', plot=False):
    """ 
    Function to calcualte RSAM of trace within a window, in minutes.
    If plot = 1, then a plot will be launched. 
    
    Parameters
    ----------
    trace: trace
        ObsPy trace to calculate RSAM on     
    window : int
        Window length in minutes.
    timePos: string
        Specifies where in the RSAM window the time is marked.  Options are
        'begin', 'middle' and 'end'. Default is middle
    plot : boolean
        Variable to determine whether or not to generate plot with results (Optional).  
        Defaults to False.
        
    Output
    ---------
    tvec: list
        list of UTCDateTimes corresponding to RSAM values.  
    rsam: list
        list of RSAM values
    
    """
    from matplotlib.dates import date2num
    import numpy as np
    t = trace.stats.starttime
    tvec = []
    rsam = []
    plotvec = []
    while t < trace.stats.endtime:
        trace_tmp = trace.slice(t, t + window * 60)
        value = np.sqrt(np.mean(np.square(trace_tmp.data)))
        if timePos == 'begin':
            stampTime = t
        elif timePos == 'middle':
            stampTime = t + (window / 2) * 60
        else:
            stampTime = t + window * 60
        plotvec.append(date2num(stampTime.datetime))
        tvec.append(stampTime)
        rsam.append(value)
        t += window * 60
    if plot == 1:
        import matplotlib.pyplot as plt
        plt.figure(figsize=[12, 3])
        plt.clf()
        plt.plot_date(plotvec, rsam, fmt='-')
        tstr = '%d min RSAM, station %s' % (window, trace.stats.station)
        plt.title(tstr)
    return tvec, rsam


def plotParticleMotion(zWave, nWave, eWave):
    """
    Plots particle motion of 3 component data.
    
    Parameters
    ----------
    zWave: Obspy stream object 
        Vertical waveform of interest for plotting particle motion.  Should be
        cut and filtered already.
    nWave: Obspy stream object
        North oriented waveform of interest for plotting particle motion.
    eWave: Obspy stream object
        East oriented waveform of interest for plotting particle motion.
        
    """

    import matplotlib.pyplot as plt
    import numpy as np
    # Start plotting particle motions
    Z = zWave
    N = nWave
    E = eWave

    # Normalize
    wave.normalize(global_max=True)

    # Create Plot
    plt.figure()

    # This does color progression, but I never got it to look right
    # npts = Z[0].stats.npts
    # cm = plt.get_cmap("RdYlGn")
    # col = [cm(float(i)/(npts)) for i in xrange(npts)]

    ax = plt.subplot(2, 2, 1)
    t = np.linspace(Z[0].stats.starttime.timestamp,
                    Z[0].stats.endtime.timestamp,
                    Z[0].stats.npts)
    plt.plot(t, Z[0].data)
    plt.title('Vertical Seismogram')
    ax = plt.subplot(2, 2, 2)
    plt.plot(E[0].data, N[0].data)
    # plt.scatter(E[0].data, N[0].data,s=10, c=col, marker='o')
    plt.xlim((-1, 1))
    plt.ylim((-1, 1))
    # plt.axis('equal')
    plt.title('Map View')
    ax = plt.subplot(2, 2, 3)
    plt.plot(N[0].data, Z[0].data)
    plt.xlim((-1, 1))
    plt.ylim((-1, 1))
    # plt.axis('equal')
    plt.title('N-S X-Section')
    ax = plt.subplot(2, 2, 4)
    plt.plot(E[0].data, Z[0].data)
    plt.xlim((-1, 1))
    plt.ylim((-1, 1))
    # plt.axis('equal')
    plt.title('E-W X-Section')

    # Do polarization
    # leigenv1, leigenv2, leigenv3, rect, plan, dleigenv, drect, dplan = eigval(E[0].data, N[0].data, Z[0].data, [1,1,1,1,1], normf=1)
    # print "Smallest eigenvalue: %f\n" % leigenv1
    # print "Intermediate eigenvalue: %f\n"  % leigenv2
    # print "Largest eigenvalue: %f\n" % leigenv3
    # print "Rectilinearity: %s\n" % rect


def streamCrossCorrelate(streamData, normalize=True, plot=True):
    """
    Takes an obspy stream object and cross-correlates the waveforms.  Typically
    good to filter before this step.
    
    Parameters
    ----------
    streamData: Obspy stream object 
        Obspy stream object, typically several events recorded on the same station.
    normalize: boolean
        If 'True', then the waveforms are normalized prior to the cross-correlation.
    plot: boolean
        If plot='True', then the cross-correlation matrix is plotted.
    
    Outputs
    ----------
    corrArray: array
        Array of cross-correlation values 
    
    """
    import matplotlib.pyplot as plt
    from obspy.signal.cross_correlation import correlate, xcorr_max
    import numpy as np

    corr = []
    ccEvt = streamData.copy()
    if normalize:
        ccEvt.normalize()
    # corr = [xcorr(ccEvt[n].data,ccEvt[m].data,500) for n in range(len(ccEvt)) for m in range(len(ccEvt))]
    corrArray = np.zeros((len(ccEvt), len(ccEvt)))
    for n in range(len(ccEvt)):
        for m in range(len(ccEvt)):
            if m > n:
                if len(ccEvt[n].data) != len(ccEvt[m].data):
                    corrArray[n][m] = 0
                    corrArray[m][n] = 0
                else:
                    cc = correlate(ccEvt[n].data, ccEvt[m].data, 200)
                    shift, corr = xcorr_max(cc)
                    corrArray[n][m] = np.abs(corr)
                    corrArray[m][n] = np.abs(corr)
            elif m == n:
                corrArray[m, m] = 1
    if plot:
        f, ax = plt.subplots()
        im = ax.imshow(corrArray, interpolation="nearest")
        f.colorbar(im, label='Cross-corrlation Value')
        f.show()

    return corrArray


def s2n(wave, t, window, method="mean", t2=None):
    """
    Calculates signal to noise ratio prior to a given time and given time window.

    Parameters
    ----------
    wave: Obspy stream object
        Obspy stream object, should be a single waveform.
    t: UTCDateTime
        Time for end of noise window and start of signal window (for example the p-wave arrival)
    window: float
        Window width in seconds to calculate "noise" and "signal" within
    method: string
        Either "mean" or "median" depending on your preference for averaging
    t2: UTCDateTime
        Time for alternate start of signal window (like the s-wave time)

    Outputs
    ----------
    sn: float
        signal to noise
    """
    import numpy as np

    if t2 == None:
        t2 = t

    # Clip noise window
    noise = wave.copy()
    noise.trim(t - window, t)
    # Clip data window
    data = wave.copy()
    data.trim(t2, t2 + window)

    if method == 'mean':
        sn = np.mean(np.abs(data.data)) / np.mean(np.abs(noise.data))
    else:
        sn = np.median(np.abs(data.data)) / np.median(np.abs(noise.data))

    return sn

def _nearest_pow_2(x):
    """
    Find power of two nearest to x

    >>> _nearest_pow_2(3)
    2.0
    >>> _nearest_pow_2(15)
    16.0

    :type x: float
    :param x: Number
    :rtype: Int
    :return: Nearest power of 2 to x
    """
    a = math.pow(2, math.ceil(np.log2(x)))
    b = math.pow(2, math.floor(np.log2(x)))
    if abs(a - x) < abs(b - x):
        return a
    else:
        return b


def checkLength(S):
    import numpy as np
    sameLength = True

    length = []
    for nowtrace in S:
        length.append(len(nowtrace))
    if np.mean(length) == length[0]:
        return sameLength
    else:
        sameLength = False
        return sameLength


def trimToShortestWindow(S):
    from obspy import UTCDateTime
    t1 = UTCDateTime(1900, 1, 1)
    t2 = UTCDateTime()
    for nowtrace in S:
        if nowtrace.stats.starttime > t1:
            t1 = nowtrace.stats.starttime
        if nowtrace.stats.endtime < t2:
            t2 = nowtrace.stats.endtime

    S.trim(t1, t2)
    return S

def medianFilter(trace, windowLength):
    import numpy as np
    datamin = np.array([])
    for windata in trace.slide(window_length=windowLength, step=windowLength):
        datamin = np.append(datamin, np.median(windata.data))
    trace.data = datamin
    trace.stats.sampling_rate = 1.0 / windowLength
    trace.stats.delta = float(windowLength)
    trace.stats.npts = len(datamin)
    trace.verify()
    

def saraCalc(trace, windowLength, decimate=False, decFactor=4, envelope=True):
    '''
    Function to calculate the median amplitude of a trace similar to Taisne (2011). Data is changed in place, so the
    original data is gone.

    :param trace: obspy Trace
    :param windowLength: length of median window in seconds
    :param decimate: boolean value on whether to decimate the trace first (suggested for long traces)
    :param decFactor: integer decimation factor
    :param envelope: boolean value on whether to take the envelope
    :return:
    '''
    
    import saraUtils
    saraUtils.saraCalc(trace, windowLength, decimate=decimate, decFactor=decFactor, envelope=envelope)
    '''
    import obspy.signal.filter as filter
    import numpy as np
    import time

    # Decimate
    if decimate:
        print ('Decimating...')
        trace.decimate(decFactor)
        print ('Done decimating...')
    # Envelope
    if envelope:
        ticenv = time.perf_counter()
        data_envelope = filter.envelope(trace.data)
        trace.data = data_envelope
        tocenv = time.perf_counter()
        print ('Envelope Elapsed Time %0.2f s' % (tocenv - ticenv,))

    # Median filter
    ticmin = time.perf_counter()
    datamin = np.array([])
    for windata in trace.slide(window_length=windowLength, step=windowLength):
        datamin = np.append(datamin, np.median(windata.data))
    trace.data = datamin
    trace.stats.sampling_rate = 1.0 / windowLength
    trace.stats.delta = float(windowLength)
    trace.stats.npts = len(datamin)
    trace.verify()
    tocmin = time.perf_counter()
    print ('Median Elapsed Time %0.2f s' % (tocmin - ticmin,))
    '''
    
def plotSpectrogram(S, figsize=[8,15], minfreq=0, maxfreq=25, tIndent=0.05, bIndent=0.05, lIndent=0.1, width=0.8,
                    waveProportion=0.25, specProportion=0.75, wavelimit=None,
                    wlen=3.0, smoothFac=2.0, per_lap=0.86, dbscale=True, clip=[0.7,1.0], userCmap='jet', lims=None, colorbar=False, ax=None,
                    fixedColor=False, vmin=-175, vmax=-110):

    '''

    :param S:
    :param figsize:
    :param minfreq:
    :param maxfreq:
    :param tIndent:
    :param lIndent:
    :param width:
    :param waveProporation:
    :param specProportion:
    :param wlen:
    :param smoothFac:
    :param per_lap:
    :param dbscale:
    :param clip:
    :param userCmap:
    :return:

    figsize = [8, 15]
    minfreq = 0
    maxfreq = 25
    tIndent = 0.05  # Top Indent
    bIndent = 0.05  # Bottom Indent
    lIndent = 0.1  # Left Indent
    width = 0.8  # Width
    waveProportion = 0.25  # Proportion of station pane to make wave
    specProportion = 0.75  # Proportion of staiton pane to make spectrogram (better add up to 1)

    # Spectrogram options
    wlen = 3.0
    smoothFac = 2.0
    per_lap = 0.86
    dbscale = True
    clip = [0.7, 1.0]
    userCmap = 'jet'
    '''

    date_format = mdates.DateFormatter('%H:%M')

    S1 = S.copy()
    #S1= sortDist(S)  # Sort stream based on distance
    nstas = len(S1)

    height = 1 - tIndent - bIndent
    staWidth = float(height) / nstas
    waveWidth = waveProportion * staWidth
    specWidth = specProportion * staWidth
    waveCursor = 1 - tIndent - waveWidth
    specCursor = 1 - tIndent - waveWidth - specWidth
    allAxes = []
    if not ax:
        #fig = plt.figure(figsize=figsize)
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_position([0, 0, 1, 1])
    else:
        fig = ax.figure
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    if fixedColor:
        norm = Normalize(vmin, vmax, clip=True)
    for tr in S1:
        if allAxes:
            #ax1 = fig.add_axes([lIndent, waveCursor, width, waveWidth], sharex=allAxes[0])
            ax1 = ax.inset_axes([lIndent, waveCursor, width, waveWidth], sharex=allAxes[0])
        else:
            #ax1 = fig.add_axes([lIndent, waveCursor, width, waveWidth])
            ax1 = ax.inset_axes([lIndent, waveCursor, width, waveWidth])
        allAxes.append(ax1)
        #t = np.arange(tr.stats.npts) / tr.stats.sampling_rate
        t = tr.times('matplotlib')
        ax1.plot(t, tr.data, 'k', linewidth=0.25, zorder=100)
        ax1.set_ylabel(tr.stats.station)
        ax1.tick_params(axis='y', length=0)
        ax1.set_yticklabels('')
        if lims:
            ax1.set_xlim(lims[0].matplotlib_date, lims[1].matplotlib_date)
        else:
            ax1.set_xlim([min(t), max(t)])
        if wavelimit:
            ax1.set_ylim([-wavelimit, wavelimit])
        ax1.get_xaxis().set_visible(False)
        #ax2 = fig.add_axes([lIndent, specCursor, width, specWidth], sharex=allAxes[0])
        ax2 = ax.inset_axes([lIndent, specCursor, width, specWidth], sharex=allAxes[0])
        allAxes.append(ax2)
        # Spectrogram generation ------------------------------------------------------------
        samp_rate = float(tr.stats.sampling_rate)
        npts = len(tr.data)
        nfft = int(_nearest_pow_2(wlen * samp_rate))
        if nfft > npts:
            nfft = int(_nearest_pow_2(npts / 8.0))
        mult = int(_nearest_pow_2(smoothFac))
        mult = mult * nfft
        nlap = int(nfft * float(per_lap))
        data = tr.data - tr.data.mean()
        end = npts / samp_rate
        specgram, freq, time = mlab.specgram(tr.data, Fs=samp_rate, NFFT=nfft,
                                             pad_to=mult, noverlap=nlap)
        # db scale and remove zero/offset for amplitude
        if dbscale:
            specgram = 10 * np.log10(specgram[1:, :])
        else:
            specgram = np.sqrt(specgram[1:, :])
        freq = freq[1:]
        if not fixedColor:
            vmin, vmax = clip
            _range = float(specgram.max() - specgram.min())
            vmin = specgram.min() + vmin * _range
            vmax = specgram.min() + vmax * _range
            norm = Normalize(vmin, vmax, clip=True)
        # calculate half bin width
        halfbin_time = (time[1] - time[0]) / 2.0
        halfbin_freq = (freq[1] - freq[0]) / 2.0
        # Setup colormap
        if userCmap:
            cmap = copy(plt.get_cmap(userCmap))
        else:
            cmap = copy(plt.cm.jet)
        cmap.set_bad(alpha=0.0)
        specgram = np.flipud(specgram)
        # center bin
        minDateNum = (tr.stats.starttime + (time[0] - halfbin_time)).matplotlib_date
        maxDateNum = (tr.stats.starttime + (time[-1] + halfbin_time)).matplotlib_date
        #extent = (time[0] - halfbin_time, time[-1] + halfbin_time,
        #          freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
        extent = (minDateNum, maxDateNum, freq[0] - halfbin_freq, freq[-1] + halfbin_freq)
        ax2.imshow(specgram, interpolation="nearest", extent=extent, norm=norm, cmap=cmap)

        if colorbar:
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), orientation='vertical', label='dB')
        #alldates = []
        #for nowt in time:
        #    #alldates.append((tr.stats.starttime+(nowt+halfbin_time)).datetime)
        #    alldates.append((tr.stats.starttime+(nowt+halfbin_time)).matplotlib_date)

        #ax2.pcolormesh(alldates, freq, specgram, norm=norm, cmap=cmap, shading='nearest')

        ax2.axis('tight')
        if lims:
            ax2.set_xlim(lims[0].matplotlib_date, lims[1].matplotlib_date)
        else:
            ax2.set_xlim(ax1.get_xlim())
        ax2.grid(False)
        # -------------------------------------------------------------------------------------

        ax2.set_ylim([minfreq, maxfreq])
        ax2.get_xaxis().set_visible(False)
        waveCursor = waveCursor - specWidth - waveWidth
        specCursor = specCursor - waveWidth - specWidth

    ax2.get_xaxis().set_visible(True)
    
    ticks = ax2.get_xticks()
    oldLabels = ax2.get_xticklabels()
    newLabels = []
    for nowtick in ticks:
        newLabels.append((S[0].stats.starttime + nowtick).strftime('%H:%M'))
    ax2.set_xticklabels(newLabels)
    
    ax2.xaxis_date()
    ax2.tick_params(axis='x', labelbottom='on')
    locator = mdates.AutoDateLocator(minticks=5, maxticks=10)
    formatter = mdates.ConciseDateFormatter(locator)
    ax2.xaxis.set_major_locator(locator)
    ax2.xaxis.set_major_formatter(formatter)

    return fig, ax1, ax2fig, ax1, ax2plot

def multiDaySpectrogram(S, averageLength=3600, fftLength=60, minFreq=0.05, maxFreq=25, cmap='magma', dateLimits=None, 
                        vmin=0.2, vmax=0.7, plotAverage=True): 


    spectraG = []
    spectraF = []
    stTimes = []
    dbscale = True
    if minFreq < 0.5:
        twoPlots = True
    else:
        twoPlots = False


    # Do the fft work
    S.merge(fill_value='interpolate')
    if dateLimits:
        S.trim(starttime=dateLimits[0], endtime=dateLimits[1])
    cntr = 1
    for wst in S.slide(window_length=averageLength, step=averageLength/10, include_partial_windows=True):  # step=averagelegnth/10 --> 3 min pixel widths/6 mins cut off on end
                                                                                                           # --> 6 mins added to spectra, 6 mins subtracted from RSAM axis
        if cntr%24 == 0:
            wst[0].stats.starttime
            #print(f"Processing window: {wst[0].stats.starttime}")
        F, G = stackSpectraContinuous(wst.detrend(type='demean'), window_length=fftLength, step=fftLength/2) # reduce fftlength for more freq resolution
                                                                                                             # (more bins, smaller height)
        spectraF.append(F)
        spectraG.append(np.abs(G))
        stTimes.append(wst[0].stats.starttime)
        cntr += 1
    #print(f"Window {cntr}: Start Time = {wst[0].stats.starttime}, End Time = {wst[0].stats.endtime}")
    #print(f"Full: Start Time = {stTimes[0]}, End Time = {stTimes[-1]}")

    # Pre-process for spectrogram
    alldates = [nowtime.matplotlib_date for nowtime in stTimes] # x dimension width for pixels
    freq = spectraF[0]                                          # argument for y dimension height for pixels
    specgram = np.flipud(np.array(spectraG).T)                  # data that fills the pixels

    # db scale and remove zero/offset for amplitude
    if dbscale:
        specgram = 10 * np.log10(specgram[:, :])
    else:
        specgram = np.sqrt(specgram[:, :])
    _range = float(specgram.max() - specgram.min())
    ##vmin = specgram.min() + vmin * _range
    ##vmax = specgram.min() + vmax * _range
    
    # standardized color scaling to our preference instead of auto data range (commented out above) #
    vmin=10
    vmax=35
    ##-------------------------------------------
    
    norm = Normalize(vmin, vmax, clip=True)
    cmap = copy(plt.get_cmap(cmap))
    cmap.set_bad(alpha=0.0)

    if plotAverage:
        Sfilt = S.copy().detrend(type='demean').taper(max_percentage=None, max_length=10).filter(type='bandpass',
                                                                                                 freqmin=minFreq,
                                                                                                 freqmax=maxFreq)
        saraCalc(Sfilt[0], 60, decimate=False, envelope=True)

    # Plot data
    f = plt.figure(figsize=(15,8))
    if minFreq < 0.5:
        if plotAverage:
            gs = gridspec.GridSpec(nrows=7, ncols=1)
            ax0 = f.add_subplot(gs[0,0])
            ax0.plot(Sfilt[0].times('matplotlib'), Sfilt[0].data, '-k')
            ax0.grid()
            ax0.get_xaxis().set_visible(False)
            #ax0.xaxis.set_ticklabels([])
            ax = f.add_subplot(gs[1:4,0], sharex=ax0)
            #ax.get_xaxis().set_visible(False)
        else:
            ax = f.add_subplot(211)
    else:
        if plotAverage:
            gs = gridspec.GridSpec(nrows=4, ncols=1)
            ax0 = f.add_subplot(gs[0, 0])
            ax0.plot(Sfilt[0].times('matplotlib'), Sfilt[0].data, '-k')
            ax0.grid()
            #ax0.get_xaxis().set_visible(False)
            #ax0.xaxis.set_ticklabels([])
            ax = f.add_subplot(gs[1:, 0], sharex=ax0)
        else:
            ax = f.add_subplot(111)

    mesh = ax.pcolormesh(alldates, np.flipud(freq), specgram, norm=norm, cmap=cmap, shading='gouraud') # sets pixel dimensions
    
    ## ----------colorbar-----------------------------------------
    # assigns line above to variable 'mesh' for colorbar
    # creates a small inset axis (colorbar scale) outside the main plot
    cax = inset_axes(ax, width="2%", height="100%", loc='lower left',
                    bbox_to_anchor=(1.02, 0., 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0)

    cbar = plt.colorbar(mesh, cax=cax)
    cbar.set_label("Power (dB)")

    # cbar = plt.colorbar(mesh, ax=ax)
    # cbar.set_label("Power (dB)")

    ##------------------------------------------------------------

    locator = mdates.AutoDateLocator(minticks=5, maxticks=15)
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)

    # setting up x-axis limits based on the range of alldates
    start_time = alldates[0] # first timestamp in data
    end_time = alldates[-1]  # last timestamp in data
    end_time_minus = end_time - (6*60)
    ax.set_xlim(start_time, end_time)

    if twoPlots:
        ax.set_ylim([1,maxFreq])
        ax.set_ylabel('Frequency, Hz')
        ax.get_xaxis().set_visible(False)
        if plotAverage:
            ax2 = f.add_subplot(gs[4:, 0], sharex=ax0)
        else:
            ax2 = f.add_subplot(212, sharex=ax)
            
        mesh2 = ax2.pcolormesh(alldates, np.flipud(freq), specgram, norm=norm, cmap=cmap, shading='gouraud')
        cax = inset_axes(ax, width="2%", height="100%", loc='lower left',
                    bbox_to_anchor=(1.02, 0., 1, 1),
                    bbox_transform=ax.transAxes,
                    borderpad=0)
        cbar = plt.colorbar(mesh2, cax=cax)
        cbar.set_label("Power (dB)")

        locator = mdates.AutoDateLocator(minticks=5, maxticks=15)
        formatter = mdates.ConciseDateFormatter(locator)
        ax2.xaxis.set_major_locator(locator)
        ax2.xaxis.set_major_formatter(formatter)
        
        # spectrogram x-limits
        ax2.set_xlim(start_time, end_time)  
        # adjust RSAM plot x-limits, stops before added time
        ax0.set_xlim(start_time, end_time_minus) # ax0 for the RSAM plot
        
        ax2.set_yscale('log')
        ax2.set_ylim([minFreq,1])
        ax2.set_xlabel('Date, UTC')
        ax2.set_ylabel('Frequency, Hz')
    else:
        ax.set_ylim([minFreq,maxFreq])
        ax.set_xlabel('Date, UTC')
        ax.set_ylabel('Frequency, Hz')

    #print('continue')

    return f


def detectClip(trace):

    clip = False
    bit = 12
    if trace.stats.channel in ['EHZ', 'EHN', 'EHE', 'EH1', 'EH2', 'SHZ', 'SHN', 'SHE'] and trace.stats.network in ['UW', 'NC']:
        tr1 = trace.copy()
        tr1.detrend(type='demean')
        maxAmpCounts = max(abs(tr1.data))
        if trace.stats.station in ['RER', 'MBW', 'GPW', 'RPW', 'RVC', 'STAR'] and trace.stats.location == '01':  # 16 bit stations
            if maxAmpCounts > 32500.0 and maxAmpCounts < 33000.0:
                clip = True
            bit = 16
        else: # 12 bit stations
            if maxAmpCounts > 2000.0 and maxAmpCounts < 2500.0:
                clip = True
    return clip, bit

def applyClipMask(S):

    from scipy.signal import find_peaks
    test = True
    S1 = S.copy()
    S1.detrend(type='demean')
    for tr in S1:
        clip, bit = detectClip(tr)
        if clip:
            thresh = 2000
            if bit == 16:
                thresh = 32500
            peaks,_ = find_peaks(abs(tr.data), height=thresh)
            times = tr.times()
            t1 = times[min(peaks)]
            t2 = times[max(peaks)]
            mask = np.zeros(len(tr.data))
            mask[min(peaks):max(peaks)+1] = 1
            boolmask = mask >= 1
            datamask = ma.MaskedArray(data=tr.data, mask=boolmask)
            if test:
                import matplotlib.pyplot as plt
                plt.subplot(2,1,1)
                plt.plot(tr.times(), tr.data, '-k')
                plt.subplot(2,1,2)
                plt.plot(tr.times(), datamask, '-k')
                plt.show()
                print('Start Time: {}'.format(t1))
                print('End Time: {}'.format(t2))

            tr.data = datamask
    return S1

def detectDeadAnalog(S):
    import numpy as np
    S1 = S.copy()
    for tr in S1:
        if tr.stats.channel not in ['BHZ', 'BHN', 'BHE', 'HHZ', 'HHN', 'HHE']:
            datamin = np.min(tr.data)
            datamax = np.max(tr.data)
            if datamax-datamin <= 50 and np.std(tr.data)<=5:
                print('{} likely dead, removing'.format(tr.id))
                S1.remove(tr)
    return S1

def getRepoData(staList, t1, t2, client=None):
    '''
    Function to get repository data, in particular from different repositories (IRIS vs. NCEDC)
    :param staList: list of NN.SSSSS.LL.CCC strings
    :param t1: UTCDateTime of the start of the window
    :param t2: UTCDateTime of the end of the window
    :param client: client to use to get data if not IRIS or NCEDC
    :return: returns an obspy stream with data and response information attached
    '''

    st = Stream()
    bulkNC = []
    bulkIRIS = []
    for nowsta in staList:
        net, sta, loc, chan = nowsta.split('.')
        if net in ['BK', 'NC', 'GM']:
            bulkNC.append((net, sta, loc, chan, t1, t2))
        else:
            bulkIRIS.append((net, sta, loc, chan, t1, t2))

    if bulkNC:
        client = Client('https://service.ncedc.org', timeout=300)
        try:
            st += client.get_waveforms_bulk(bulkNC, attach_response=True)
        except Exception as error:
            print('No data from NCEDC, or other error.')
            print("The actual error:", error)
    if bulkIRIS:
        client = Client('IRIS', timeout=300)
        try:
            st += client.get_waveforms_bulk(bulkIRIS, attach_response=True)
        except Exception as error:
            print('No data from IRIS, or other error.')
            print("The actual error:", error)

    return st

def preAmpProcessing(staList, t1, t2, staFile, saraint=30, snr=None, dataFile=None, removeClip=True, checkDead=True, 
                     correctionFile=None, noiseFile=None, pre_filt=(0.5,1,10,11), realTime=False, config=None, returnRaw=False):
    
    import amp
    import ampCorrect
    import intensityLocate
    # Get relative times to account for noise windows and signal processing
    noiseStart = t1 - 600
    noiseEnd = noiseStart + 600

    dataStart = t1 - 60
    dataEnd = t2 + 60

    # ---------------------------------------------------------------
    # Get station information
    chanList, chans = amp.readStaInfo(staFile)
    if dataFile:
        S = read(dataFile)  # Datafile should have bad data removed
    else:
        st0 = Stream()
        if not realTime:
            st = getRepoData(staList, noiseStart, dataEnd)
        else:
            st = getRepoData(staList, dataStart, dataEnd)
        if not st:
            print('No data')
            return Stream(), [], chans
        if removeClip: # Check for clipping
            for tr in st:
                clip, bit = detectClip(tr)
                # clip=False
                if not clip:
                    st0.append(tr)
                else:
                    print('{} discarded because clip detected!'.format(tr.id))
        else:
            st0 = st.copy()
        
        if checkDead: # Detect dead analog
            st00 = detectDeadAnalog(st0)
            st1 = st00.copy()
        st1.remove_response(output='VEL', pre_filt=pre_filt)
        st2 = st1.copy()
        for tr in st2:
            if tr.stats.sampling_rate > 51:
                tr.decimate(factor=2)
        st2.merge(fill_value=0)
        if not realTime:
            st2.trim(starttime=dataStart, endtime=dataEnd, pad=True, fill_value=0)
        else:
            st2.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
        S = st2.copy()
    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    # Get channels in appropriate order
    S1 = Stream()
    for nowChan in chanList:
        net, sta, loc, chan = nowChan.split('.')
        S1 = S1 + S.select(network=net, station=sta, location=loc, channel=chan)
    waveRaw = S1.copy()
    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    # Create envelope    
    for tr in S1:
        saraCalc(tr, saraint)
    envRaw = S1.copy()
    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    # Signal to Noise cutoff
    if snr:
        if realTime:
            if noiseFile:
                print('noiseFile defined in two places!  Going with the config file!')

            S2, SNtrace, sn = intensityLocate.applySignalToNoiseMaskRT(S1, config)
        elif noiseFile:
            path = Path(noiseFile)
            if path.is_file():
                print('Using existing noiseFile: {}'.format(noiseFile))
                S2, SNtrace, sn = intensityLocate.applySignalToNoiseMaskFile(S1, noiseFile, snRatio=snr)
            else:
                S2, SNtrace, sn = intensityLocate.applySignalToNoiseMask(S1, noiseStart, noiseEnd, snRatio=snr, writeFile=noiseFile)
        else:
            S2, SNtrace, sn = intensityLocate.applySignalToNoiseMask(S1, noiseStart, noiseEnd, snRatio=snr)
    else:
        S2 = S1.copy()
        sn = {}
        SNtrace = S1.copy()
        for tr in SNtrace:
            tr.data = tr.data*0
    # ---------------------------------------------------------------
    # ---------------------------------------------------------------
    # Apply amplitude correction
    if correctionFile:
        S2, corrs = ampCorrect.applyCorrectionFile(S2, correctionFile)
    else:
        print('No correction file given, no corrections applied.')

    # ---------------------------------------------------------------
    # Get times to use
    tAll = S1[0].times("utcdatetime")  # Use the stream before it is masked
    t = []
    for tNow in tAll:
        if tNow > t1 and tNow < t2:
            t.append(tNow)

    if returnRaw:
        if correctionFile:
            corrwaves,_ = ampCorrect.applyCorrectionFile(waveRaw, correctionFile)
            correnv,_ = ampCorrect.applyCorrectionFile(envRaw, correctionFile)
        else:
            corrwaves = waveRaw
            correnv = envRaw
        return S2, t, chans, corrwaves, correnv
    else:
        return S2, t, chans