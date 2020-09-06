#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Lime Auto Fast
# GNU Radio version: 3.7.14.0
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from baz import facsink
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import wxgui
from gnuradio.eng_option import eng_option
from gnuradio.fft import window
from gnuradio.filter import firdes
from gnuradio.wxgui import fftsink2
from grc_gnuradio import wxgui as grc_wxgui
from optparse import OptionParser
import limesdr
import wx


class Lime_Auto_Fast(grc_wxgui.top_block_gui):

    def __init__(self):
        grc_wxgui.top_block_gui.__init__(self, title="Lime Auto Fast")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 50e3

        ##################################################
        # Blocks
        ##################################################
        self.wxgui_fftsink2_0 = fftsink2.fft_sink_c(
        	self.GetWin(),
        	baseband_freq=0,
        	y_per_div=10,
        	y_divs=10,
        	ref_level=0,
        	ref_scale=2.0,
        	sample_rate=samp_rate,
        	fft_size=1024,
        	fft_rate=15,
        	average=False,
        	avg_alpha=None,
        	title='FFT Plot',
        	peak_hold=False,
        )
        self.Add(self.wxgui_fftsink2_0.win)
        self.limesdr_source_0 = limesdr.source('', 0, '')
        self.limesdr_source_0.set_sample_rate(samp_rate)
        self.limesdr_source_0.set_center_freq(393.5e6, 0)
        self.limesdr_source_0.set_bandwidth(5e6,0)
        self.limesdr_source_0.set_gain(39,0)
        self.limesdr_source_0.set_antenna(3,0)
        self.limesdr_source_0.calibrate(5e6, 0)

        self.facsink_0 = facsink.fac_sink_c(
        	self.GetWin(),
        	title='Fast AutoCorrelation',
        	sample_rate=samp_rate,
        	baseband_freq=0,
                y_per_div=5,
        	ref_level=-40,
        	fac_size=65536,
                fac_rate=facsink.default_fac_rate,
                average=True,
        	avg_alpha=0.3,
        	peak_hold=True,
        )
        self.Add(self.facsink_0.win)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.limesdr_source_0, 0), (self.facsink_0, 0))
        self.connect((self.limesdr_source_0, 0), (self.wxgui_fftsink2_0, 0))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.wxgui_fftsink2_0.set_sample_rate(self.samp_rate)
        self.facsink_0.set_sample_rate(self.samp_rate)


def main(top_block_cls=Lime_Auto_Fast, options=None):

    tb = top_block_cls()
    tb.Start(True)
    tb.Wait()


if __name__ == '__main__':
    main()
