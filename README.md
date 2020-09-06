# LimeAuto

Introduction

Lime-Auto is a real time implementation of the AutoCorrelation Function (AFC) on RF signal for LimeSDR radios. It allows Automatic Signal Identification (ASI) in a few seconds whatever is the modulation used.

LimeAuto is a modified version of LimeScan belonging to Lime-Tools suite (https://github.com/myriadrf/lime-tools) to implement Automatic Signal Identification (ASI) capacity.

The current implementation is working for narrowband signal (TETRA, TETRAPOL, MPT1237, POCSAG, DMR, INMARSAT Standard C, and STANAG4285) as well as large band signal (GSM, LTE and DAB).

Analysis is either on a specific given frequency or within a specified frequency band.

See more details on my YouTube channel (sdrpower2) with two videos on ACF (https://www.youtube.com/watch?v=5Oc8XUg7aFo&t=28s) and ASI (.
https://www.youtube.com/watch?v=O9Qmykk6N4U&t=24s)

Installation

cd lime-auto
mkdir build
cd build
cmake ..
make
sudo make install
sudo ldconfig

Command line

To do identification through on a specific frequency (here 466.050M)
sudo LimeAuto -f 466.050M:466.050M -ns 1 -d 0 -g 29

To do identification through a frequency band (here 120 * 10 kHz = 1,2 MHz band  starting at 393.250 MHz)
sudo LimeAuto -f 393.250M:394.450M -ns 120 -d 0 -g 29 -st 10k

With the additional parameter -verb (0 or 1 or 2) you can have more or lesss details on the signal identification parameters outputs.

Key and expert parameters

The first key parameter is the sample rate  (Fs= 50 kHz). This parameter is not to be changed as all the algoritms are based on this value and will not work for another one.  As a consequence, the modification of this parameter from the command line has been disabled.

The second key parameter is the FFT length (NFFT= 32768 is the recommended), which is a good compromise between speed and time depth of the AutoCorrelation Function. As it is recommended not to change it, the modification of this parameter from the command line has been disabled.

Keep in mind the following relationship about ACF time depth

			ACF time depth = (NFFT/2)  / Fs

For the above parameters (50 kHz and 32 768), ACF time length = 327 ms. This means that if the recurrence of the autocorrelation peak is 106.66 ms (STANAG 4285), the maximum number of peaks that you can detect is 3. You can increase it by increasing NFFT but processing time will increase too.

The third key parameter is the gain, which has not to be too low or to be too high. This is obviously depending on the environment (spectral environment, antenna used, LNA or not). The modification of this parameter is accessible from the command line (- g).

The last key parameter is relative to the information given about the results of the identification algorithm (-verb). You can choose between 0 (just the result), 1 (which gives you some idea of the ranking between all the signals candidate) and 2 (more details).

The main other parameters are “expert parameters”. They can be only be modified through recompilation. Just a few words about it:
•	The number of FFT for RMS integration. 5 is a good compromise between processing speed and to have some RMS integration
•	The mutiplying factor used with the standard deviation to filter what is considered as peak or not. 1.2 is a good compromise between to filter peaks which are noise and peaks which are signals
•	The parameter to consider a  signal detection. 51 (over 100) is considered as a good compromise between false detection and missed detections
•	The number of analysis for a second RMS integration level. As there is an infinite loop until detection, there is automatically a second RMS integration level and this parameter is let to 1.

Other signals

In case no signal is recognized, it is possible to launch under gnuradio the "autocorrelation_lime_gui.grc" file, which is an ACF implementation. In case of a framed signal, autocorrelation peaks will appear and the recurrency can be measured visually.

Licence

This software is published under the Apache License 2.0.

