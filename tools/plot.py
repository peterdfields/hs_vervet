import matplotlib
import numpy
import pylab

#plot sfs (from dadi)

class _sfsTickLocator(matplotlib.ticker.Locator):
    def __call__(self):
        'Return the locations of the ticks'

        try:
            vmin, vmax = self.axis.get_view_interval()
            dmin, dmax = self.axis.get_data_interval()
        except AttributeError:
            self.verify_intervals()
            vmin, vmax = self.viewInterval.get_bounds()
            dmin, dmax = self.dataInterval.get_bounds()

        tmin = max(vmin, dmin)
        tmax = min(vmax, dmax)

        return numpy.array([round(tmin)+0.5, round(tmax)-0.5])


#: Custom tick formatter
_ctf = matplotlib.ticker.FuncFormatter(lambda x,pos: '%i' % (x-0.4))
class Sfs(numpy.ma.masked_array):
    def __new__(subtype, data, mask=numpy.ma.nomask, mask_corners=True,
                data_folded=None, check_folding=True, dtype=float, copy=True,
                fill_value=numpy.nan, keep_mask=True, shrink=True,
                pop_ids=None, extrap_x=None):
        data = numpy.asanyarray(data)

        if mask is numpy.ma.nomask:
            mask = numpy.ma.make_mask_none(data.shape)

        subarr = numpy.ma.masked_array(data, mask=mask, dtype=dtype, copy=copy,
                                       fill_value=numpy.nan, keep_mask=True,
                                       shrink=True)
        subarr = subarr.view(subtype)
        if pop_ids is not None:
            subarr.pop_ids = pop_ids
        else:
            subarr.pop_ids = ["1","1"]
        return subarr
def plot_single_2d_sfs(sfs, vmin=None, vmax=None, ax=None,
                       pop_ids=None, extend='neither', colorbar=True):
    """
    Heatmap of single 2d SFS..
....
    If vmax is greater than a factor of 10, plot on log scale.

    sfs: SFS to plot
    vmin: Values in sfs below vmin are masked in plot.
    vmax: Values in sfs above vmax saturate the color spectrum.
    ax: Axes object to plot into. If None, the result of pylab.gca() is used.
    pop_ids: If not None, override pop_ids stored in Spectrum.
    extend: Whether the colorbar should have 'extension' arrows. See
            help(pylab.colorbar) for more details.
    colorbar: Should we plot a colorbar?
    """
    if ax is None:
        ax = pylab.gca()

    if vmin is None:
        vmin = sfs.min()
    if vmax is None:
        vmax = sfs.max()

    pylab.cm.hsv.set_under('w')
    if vmax / vmin > 10:
        # Under matplotlib 1.0.1, default LogFormatter omits some tick lines.
        # This works more consistently.
        norm = matplotlib.colors.LogNorm(vmin=vmin*(1-1e-3), vmax=vmax*(1+1e-3))
        format = matplotlib.ticker.LogFormatterMathtext()
    else:
        norm = matplotlib.colors.Normalize(vmin=vmin*(1-1e-3),
                                           vmax=vmax*(1+1e-3))
        format = None
    mappable=ax.pcolor(numpy.ma.masked_where(sfs<vmin, sfs),
                       cmap=pylab.cm.hsv, shading='flat',
                       norm=norm)
    #hs
    try:
        ax.figure.colorbar(mappable, extend=extend, format=format)
    except:
        pass
    #/hs
    #original:
    #ax.figure.colorbar(mappable, extend=extend, format=format)
    if not colorbar:
        del ax.figure.axes[-1]

    ax.plot([0,sfs.shape[1]],[0, sfs.shape[0]], '-k', lw=0.2)

    if pop_ids is None:
        if sfs.pop_ids is not None:
            pop_ids = sfs.pop_ids
        else:
            pop_ids = ['pop0','pop1']
    #hs:
    ax.set_ylabel(pop_ids[0], horizontalalignment='right')
    ax.yaxis.set_label_coords(-0.05, 0.6)
    #/hs
    #ax.set_ylabel(pop_ids[0], horizontalalignment='left')
    ax.set_xlabel(pop_ids[1], verticalalignment='bottom')

    ax.xaxis.set_major_formatter(_ctf)
    ax.xaxis.set_major_locator(_sfsTickLocator())
    ax.yaxis.set_major_formatter(_ctf)
    ax.yaxis.set_major_locator(_sfsTickLocator())
    for tick in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
        tick.set_visible(False)

    ax.set_xlim(0, sfs.shape[1])
    ax.set_ylim(0, sfs.shape[0])



