/*
 Copyright 2018 Lime Microsystems Lt

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*LimeAuto is a modified version of LimeScan belonging 
to Lime-Tools suite (https://github.com/myriadrf/lime-tools)
to implement automatic signal identification (ASI) 
based on Auto Correlation Function (ACF).
Same licences conditions as for LimeScan apply.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h> // double _Complex, compliant with double, fftw_complex[2] Re=0, Im=1
#include <fftw3.h>  // should come after complex.h
#include <unistd.h> // pipes usleep
#include <fcntl.h> // file control
#include <sys/types.h> // pipes
#include <sys/stat.h> // mkdir
#include <time.h> // tools for telling the time!
#include <errno.h> // EPERM etc
#include "LimeSuite.h"

#define PI 3.14159265
#define PI2 (2*PI)


using namespace std;

void DecCmdLine( int argc, char *argv[] );
void Init( void );
void OpenSDR( void );
void Scan( void );
void GnuPlotDispAndSave( char *fName,double **zz );
void DisplayBinFile( char *fname );
void CloseSDR( void );
int error( void );
void ShutDown( void );
int ReadPwl( char fname[],double xlim[],double mPwl[],double cPwl[] );
void ApplyPwl( double *vecy[],double xlim[],double mPwl[],double cPwl[],unsigned char nPwl );
void ReadNMEAtraffic( char *gpsPath );
char* ReadTil( char *buffer, char word[] );
void ReadSubWord(char word[],char subWord[],int pos,int len);
void ReadGPRMC( char *buffer );

// default global control variables
float frq_min=800.0e6; // 800 MHz
float frq_max=1000.0e6; // 1000 MHz
float frq_step=5e4; // 30 MHz avec 192 (75%)  - 32,5 Mhz avec 208 (81,25%)
float scan_step=10e3; //
unsigned int NFFTcro=8192*4;
double frq=0.0;
float frq_LPF=5e6; // IF bandwidth of LimeSDR
float frq_Samp=5e4; // USB Sample Rate for LimeSDR
float frq_Tic=10.0;
unsigned char OSR=0; // ADC runs faster than USB rate.
unsigned int gaindB=29; // gain in dB 0:73 - int for library compatibility
unsigned char Ch=0; // SDR Channel number. 0:1
unsigned int NFFT=8192*4; // default number of FFT bins
unsigned int NRpt=5; // number of repeat reads for RMS integration
unsigned int NScan=400; // number of scan
unsigned int Dev=0; // device number
unsigned int ScanSpeed=0; // Scan Speed Measurement
unsigned int StopGo=0; // Stop Go Measurement
unsigned int Graph=0; // Graph (brut)= 1 - Graph (moyenne = 2) - Graph (pic=3)
unsigned int Verbose=1; // Verbose = 0 (no) - Verbose = 1 (low) - Verbose = 2 (moderate)
int ch;
time_t t0;
time_t t1;

// choix de l'antenne
char LNAnum=3; // off - use default LNAW = 3 LNAH = 1
// choix de l'antenne

// private global variables
char fNameStem[50]="output"; // default output location
unsigned int NUSB=1420*10; // number of samples read per USB packet
unsigned int Nlat=NUSB/NFFT+((NUSB%NFFT)>0); // ceil(NUSB/NFFT)
unsigned int fCnt;
lms_device_t* device=NULL; // SDR device
lms_stream_t streamId; // SDR stream
double *win=NULL; // hamming window coefficients
double _Complex *in=NULL; // registered with FFTW
double _Complex *out=NULL; // registered with FFTW

//double *ini=NULL; // registered with FFTW INVERSE
//double _Complex *outi=NULL; // registered with FFTWi

double *ini=NULL;
double _Complex *outi=NULL;

fftw_plan pfft; // FFTW

fftw_plan pffti; // FFTW INVERSE INVERSE

time_t *swpTime;
double **fmat; // freq values for interpolation with pwl (MHz)
double **zrms; // average data
double **zmin; // min data
double **zmax; // max data
double RNFFTdB;
char gpsPath[20]="/dev/ttyACM0";
unsigned char doGPS=0; // 0=don't 1=do
unsigned int tstLvl=-50; // for optional test signal
double tstFrq=860e6; // 860MHz is an amature band, dont want to jam sensitive bands
unsigned char doTst=0; // 0=don't 1=do, test signal defaul off.

char doPwlLNAW=0;
char doPwlLNAH=0;
char doPwlLNAL=0;
char doPwlAnt=0;
char fPwlLNAW[70]="pwl/lnaw.txt";
char fPwlLNAH[70]="pwl/lnah.txt";
char fPwlLNAL[70]="pwl/lnal.txt";
char fPwlAnt[70]="pwl/ant.txt";

int main( int argc, char *argv[] )
{
	unsigned int ct;
	unsigned char nPwl;
	double xPwl[255];
	double mPwl[255];
	double cPwl[255];
	char fname[100];
	DecCmdLine(argc,argv);
	if( doGPS>0 )
		ReadNMEAtraffic( gpsPath );
	Init();
	OpenSDR();
	if (ScanSpeed != 0)
		t0 = time(NULL); //permet de mesurer le temps donc la vitesse de scan mais faut lancer LimeScan uniquement
	frq=frq_min-scan_step;
	for (ct=0;ct<NScan;ct++)
	{
	frq=frq+(scan_step);
	//printf("\033[34;01mFréquence: %.3f\033[00m\n",frq);
	if (Graph == 0)
		printf("Frequency : %i\n",ct+1);
	Scan();

	if( (doPwlLNAH>0) && (LNAnum==1) ) // Apply PWL corrections
	{ // block pwl if for wrong LNA
		nPwl=ReadPwl( fPwlLNAH,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( (doPwlLNAL>0) && (LNAnum==2) )
	{
		nPwl=ReadPwl( fPwlLNAL,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( (doPwlLNAW>0) && (LNAnum==3) )
	{
		nPwl=ReadPwl( fPwlLNAW,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	}
	if( (doPwlAnt>0) )
	{
		nPwl=ReadPwl( fPwlAnt,xPwl,mPwl,cPwl );
		ApplyPwl( zrms,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmin,xPwl,mPwl,cPwl,nPwl );
		ApplyPwl( zmax,xPwl,mPwl,cPwl,nPwl );
	} // display results	sprintf(fname,"%sRMS",fNameStem);
//	sprintf(fname,"%sRMS",fNameStem);
	GnuPlotDispAndSave( fname,zrms ); // average
//	sprintf(fname,"%sMin",fNameStem);
//	GnuPlotDispAndSave( fname,zmin ); // max
//	sprintf(fname,"%sPk",fNameStem);
//	GnuPlotDispAndSave( fname,zmax ); // min
#ifdef USE_GNUPLOT
//	sprintf(fname,"%s/%sRMS.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname ); // provides updated window
//	sprintf(fname,"%s/%sPk.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname );
//	sprintf(fname,"%s/%sMin.gif",fNameStem,fNameStem);
//	DisplayBinFile( fname );
#endif

	}
	

	if (ScanSpeed!=0)
		{
		printf("Scan duration: %lu seconds\n",(time(NULL)-t0)); //permet de mesurer le temps donc la vitesse de scan 
		printf("Scan speed: %.0f MHz/s\n",frq_step*fCnt*NScan*1.0e-6/(time(NULL)-t0));
		}
	CloseSDR();
	ShutDown();
	return(0);
}

void Init( void )
{ // memory allocation and frequency settings
	unsigned int ci;
	in=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	out=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT); // align for SIMD
	win=(double *)malloc(sizeof(double)*NFFT);

	//ini=(double *)fftw_malloc(sizeof(double)*NFFT);
	//outi=(double _Complex *) fftw_malloc(sizeof(double _Complex)*NFFT);

	ini=(double *)fftw_malloc(sizeof(double)*NFFT);
	outi=(double _Complex*) fftw_malloc(sizeof(double _Complex)*NFFT);

	//frq_step=frq_Samp;
//printf("%.0f\n",frq_min);
//printf("%.0f\n",frq_step*0.5);

//frq_min=frq_min+frq_step/2;



	frq_max-=frq_step/2;
//printf("%.0f\n",frq_min);
	
//fCnt=ceil((frq_max-frq_min)/frq_step)+1;

//ligne à décommenter en mode tableau ou Graphique
	if (Graph != 0)	
		fCnt=1;
	if (Graph == 0)	
		fCnt=1;
//ligne à décommenter en mode tableau ou Graphique

	printf("fCnt=%i\n",fCnt);
	printf("NFFT=%i\n",NFFT);
	printf("frequence_sampling=%.0f\n",frq_Samp);
	swpTime=(time_t*)malloc(sizeof(time_t)*fCnt); // was NFFT
	fmat=(double**)malloc(sizeof(double*)*fCnt);
	zrms=(double**)malloc(sizeof(double*)*fCnt);
	zmin=(double**)malloc(sizeof(double*)*fCnt);
	zmax=(double**)malloc(sizeof(double*)*fCnt);
	for(ci=0;ci<fCnt;ci++)
	{
		fmat[ci]=(double*)malloc(sizeof(double)*NFFT);
		zrms[ci]=(double*)malloc(sizeof(double)*NFFT);
		zmin[ci]=(double*)malloc(sizeof(double)*NFFT);
		zmax[ci]=(double*)malloc(sizeof(double)*NFFT);
	}
	for(ci=0;ci<NFFT;ci++) // calcualte Windowing function coefficients
		win[ci]=1-cos((PI2*ci)/(NFFT-1)); // Hann
//		win[ci]=25.0/46.0-(1-25.0/46.0)*cos((PI2*ci)/(NFFT-1)); // Hamming
//		win[ci]=0.355768-0.487396*cos((PI2*ci)/(NFFT-1))+0.144232*cos((2*PI2*ci)/(NFFT-1))-0.012604*cos((3*PI2*ci)/(NFFT-1)); // Nutall
//		win[ci]=0.3635819-0.4891775*cos((PI2*ci)/(NFFT-1))+0.1365995*cos((2*PI2*ci)/(NFFT-1))-0.0106411*cos((3*PI2*ci)/(NFFT-1)); // Blackman Nutall
//		win[ci]=0.35875-0.48829*cos((PI2*ci)/(NFFT-1))+0.14128*cos((2*PI2*ci)/(NFFT-1))-0.01168*cos((3*PI2*ci)/(NFFT-1)); // Blackman Harris
	pfft=fftw_plan_dft_1d(NFFT,
                reinterpret_cast<fftw_complex*>(in),
                reinterpret_cast<fftw_complex*>(out),
                FFTW_FORWARD,
                FFTW_ESTIMATE);
//	pffti=fftw_plan_dft_r2c_1d(NFFT,
//                reinterpret_cast<double*>(ini),
//                reinterpret_cast<double _Complex*>(outi),
//                FFTW_ESTIMATE);
	pffti=fftw_plan_dft_r2c_1d(NFFT,
                reinterpret_cast<double*>(ini),
                reinterpret_cast<fftw_complex*>(outi),
                FFTW_ESTIMATE);
}

void ShutDown( void )
{ // deallocate memory, and other good house keeping
	unsigned int ci;
	for(ci=0;ci<fCnt;ci++)
	{
		free(fmat[ci]);
		free(zrms[ci]);
		free(zmin[ci]);
		free(zmax[ci]);
	}
	free(win);
	free(swpTime);
	fftw_destroy_plan(pfft);
	fftw_free(in);
	fftw_free(out);

	fftw_destroy_plan(pffti);
	fftw_free(ini);
	fftw_free(outi);
}

void OpenSDR( void )
{ // software based SDR set up, allows use with LimeSDR and LimeSDRmini.
	int n; // Find devices
	int ci=0;
	unsigned int gaindBTx;
	float_type rate=frq_Samp;
	float_type rf_rate=frq_Samp*OSR;
	float_type gain;

lms_name_t antenna_list[10]; // large enough list for antenna names.
	if((n=LMS_GetDeviceList(NULL))<0) // Pass NULL to only obtain number of devices
		error();
	printf("Devices found: %i \n",n);
	if(n<1)
		error();
	lms_info_str_t* list = new lms_info_str_t[n]; // allocate device list storage
	if(LMS_GetDeviceList(list)<0)
		error();
	for(ci=0;ci<n;ci++) // print device list
		printf("%i:%s\n",ci,list[ci]);
	if(LMS_Open(&device,list[Dev],NULL)) //Open device Dev
		error();
	delete [] list;
    // Initialize device with default configuration
    // Use LMS_LoadConfig(device, "/path/to/file.ini") to load config from INI
    if(LMS_Init(device) != 0)
        error();
    if(LMS_EnableChannel(device, LMS_CH_RX,Ch,true)!=0) // Rx, Channels 0:(N-1)
        error();


     if(LMS_SetLOFrequency(device, LMS_CH_RX,Ch,frq_min)!=0)
        error();
	// Alternatively, NULL can be passed to LMS_GetAntennaList() to obtain antennae num
	if((n = LMS_GetAntennaList(device, LMS_CH_RX, 0, antenna_list)) < 0)
		error();
	printf("Ae:\n"); // print available antennae names
	for(ci = 0; ci < n; ci++)
		printf(" %i:%s",ci,antenna_list[ci]);
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) // get selected antenna index
		error();
	printf("\nDefault: %i:%s, ",n,antenna_list[n]);
	if( LNAnum==1 )
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAH) != 0)
			error();
	if( LNAnum==2 )
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAL) != 0)
			error();
	if( LNAnum==3 )
		if(LMS_SetAntenna(device, LMS_CH_RX, Ch, LMS_PATH_LNAW) != 0)
			error();
	if((n = LMS_GetAntenna(device, LMS_CH_RX, Ch)) < 0) // get selected antenna index
		error();
	printf("Selected: %i:%s\n",n,antenna_list[n]);

    if(LMS_SetSampleRate(device,frq_Samp,OSR) != 0) // oversampling in RF OSR x sample rate
        error();

if(LMS_GetSampleRate(device, LMS_CH_RX, 0, &rate, &rf_rate) != 0)  // can pass NULL
		error();
	printf("\nUSB rate: %f MS/s, ADC rate: %fMS/s\n", rate*1.0e-6, rf_rate*1.0e-6);
    //Example of getting allowed parameter value range
    //There are also functions to get other parameter ranges (check LimeSuite.h)
	lms_range_t range; //Get allowed LPF bandwidth range
	if(LMS_GetLPFBWRange(device,LMS_CH_RX,&range)!=0)
		error(); // RX LPF range: 1.400100 - 130.000000 MHz
//	printf("RX LPF bandwitdh range: %f - %f MHz\n", range.min*1.0e-6,range.max*1.0e-6);
	if(LMS_SetLPFBW(device, LMS_CH_RX, Ch, frq_Samp)!=0)  //Configure LPF, bandwidth 8 MHz
		error();
	if(LMS_SetNormalizedGain(device,LMS_CH_RX,Ch,gaindB/73.0)!=0) //Set RX gain
		error();
	int err=0;
	if(LMS_SetGaindB(device,LMS_CH_RX,Ch,gaindB)!=0) // 0:73
	{
		printf("Ch=%i gaindB=%i err=%i\n",Ch,gaindB,err);
		error();
	}
	if(LMS_GetNormalizedGain(device,LMS_CH_RX,Ch,&gain)!=0) //normalized gain
		error();
	if(LMS_GetGaindB(device,LMS_CH_RX,Ch,&gaindB)!=0)
		error();
	printf("Normalized RX Gain: %f, RX Gain: %i dB\n",gain,gaindB);


	if (frq_step < 2.5e6)
	{
		//printf("Warning: Calibrating for 2.5 MHz bandwidth (requested %.2f MHz [out of range])\n", frq_step/1e6);
		if(LMS_Calibrate(device,LMS_CH_RX,Ch,2.5e6,0)!=0)
			error();
	}
	else
	{
		if(LMS_Calibrate(device,LMS_CH_RX,Ch,frq_step,0)!=0)
			error();
	}

	if( doTst>0 )
	{ // optional TX test signal
		if(LMS_EnableChannel(device,LMS_CH_TX,Ch,true)!=0)
			error();
		if(LMS_SetAntenna(device,LMS_CH_TX,Ch,LMS_PATH_TX1)!=0)
			error();
		if(LMS_SetLOFrequency(device,LMS_CH_TX,Ch,tstFrq)!=0)
			error();
		if(LMS_SetTestSignal(device,LMS_CH_TX,Ch,LMS_TESTSIG_NCODIV8,0,0)!=0)
			error(); // freq offset = rate/NCODIV
		if(LMS_SetLPFBW(device,LMS_CH_TX,Ch,16.0E6)!=0)
			error(); // TX LPF range: 1.400100 - 130.000000 MHz
		if(LMS_SetGaindB(device,LMS_CH_TX,Ch,(79+tstLvl))!= 0) // 0:73
			error(); //tstLvl
		if(LMS_GetNormalizedGain(device,LMS_CH_TX,Ch,&gain)!=0) //normalized gain
			error();
		if(LMS_GetGaindB(device,LMS_CH_TX,Ch,&gaindBTx)!=0)
			error();
		printf("Normalized TX Gain: %f, TX Gain: %i dB\n",gain,gaindBTx);
		if(LMS_Calibrate(device,LMS_CH_TX,Ch,frq_step,0)!=0)
			error();
	}
	streamId.channel=Ch; //channel number
	streamId.fifoSize=1024*1024; //fifo size in samples
	streamId.throughputVsLatency=1.0; //optimize for max throughput
	streamId.isTx=false; //RX channel
	streamId.dataFmt=lms_stream_t::LMS_FMT_F32; //32-bit floats
	if(LMS_SetupStream(device,&streamId)!=0)
		error();
	LMS_StartStream(&streamId);
}

void CloseSDR( void )
{
	if(LMS_EnableChannel(device,LMS_CH_TX,Ch,false)!=0)
		error();
	if(LMS_SetGaindB(device,LMS_CH_TX,Ch,0)!= 0) // 0:73
		error(); // switch off Tx so not to interfere.
	LMS_StopStream(&streamId); // stream is stopped, start again with LMS_StartStream()
	LMS_DestroyStream(device, &streamId); //stream can no longer be used
	LMS_Close(device);
}

int error( void )
{
	printf("LimeSDR: ERROR\n");
	if(device!=NULL)
		LMS_Close(device);
	exit(-1);
}
void Scan( void )
{
	int samplesRead;
	unsigned int ci,cf,ct,ctt,ck,bord;
	float _Complex *buf=(float _Complex*)malloc(sizeof(float _Complex)*NFFT); // LimeSDR
	float *mb[4]; // 0 mag, 1 min mag, 2 max mag, 3 rms mag, do in vector form for speed
	unsigned int NFFTd2=NFFT>>1;
	float RNRptdB=-10*log10(NRpt); // divide Sum( mag^2 )/NRpt in dBs
	float RNFFTdB=-20*log10(NFFT); // divide by NFFT in dBs
	double RFstep=frq_Samp/NFFT;
	float ecart;
	float qecart;
	int TETRAPOL_C,TETRA_C,MPT1237_C,DMR_C,POCSAG_C,STANDARDC_C,DVBT_C,GSM_C,DAB_C,LTE_C,compteur_DMR;
	float TETRAPOL_L,TETRA_L,MPT1237_L,DMR_L,POCSAG_L,STANDARDC_L,DVBT_L,GSM_L,DAB_L,seuil_DMR, LTE_L, max, winner, moyenne,somme,score;
	
	for( ci=0;ci<4;ci++)
		mb[ci]=(float *)malloc(sizeof(float)*NFFT);

	TETRAPOL_C=0, TETRAPOL_L=0;
	TETRA_C=0, TETRA_L=0;
	DMR_C=0, DMR_L=0;
	MPT1237_C=0, MPT1237_L=0;
	POCSAG_C=0, POCSAG_L=0;
	STANDARDC_C=0, STANDARDC_L=0;
	DVBT_C=0, DVBT_L=0;
	GSM_C=0, GSM_L=0;
	DAB_C=0, DAB_L=0;
	LTE_C=0, LTE_L=0;
	compteur_DMR=0; seuil_DMR=6;

	t0 = time(NULL);// to do time measurement
	
	if (StopGo == 1) 
		{
		printf("Ready?\n");
		ch = getc(stdin);
		}

// loop do-while (A) (start)  ************************************************************************************
	
	winner =0;
	do
	{

// AFC analysis (B) (start)  *************************************************************************************

	for( cf=0; cf<fCnt; cf++ )
	{

//ligne à décommenter en mode tableau ou Graphique
	if (Graph != 0)
		frq=frq_min+cf*frq_step;
//	if (Graph ==0)
//		frq=frq_min;
//		frq=frq_min+cf*frq_step;
//ligne à décommenter en mode tableau ou Graphique


		LMS_StopStream(&streamId);
		if(LMS_SetLOFrequency(device,LMS_CH_RX,0,frq)!= 0)
			error();
		//usleep(100); // flush old data, was 100us
		LMS_StartStream(&streamId);
		usleep(100); // flush old data, was 100us
		for( ci=0;ci<NFFT;ci++) // vectorized
		{ // SIMD vector reset
			mb[0][ci]=0.0; // now
			mb[1][ci]=0.0; // rms
			mb[2][ci]=0.0; // min
			mb[3][ci]=0.0; // max
		}

// First FFT (complex to complex)
		for( ct=0; ct<NRpt; ct++ )
		{
			swpTime[cf]=time(NULL);
			samplesRead=LMS_RecvStream(&streamId,buf,NFFT,NULL,NFFT);

			for(ci=0;ci<NFFT;ci++) // copy SDR buffer to FFTW and no window (Vectorized)
				in[ci]=buf[ci]; // no windowing

			fftw_execute(pfft);

			for( ci=0;ci<NFFT;ci++) // vectorized
			{
				mb[0][ci]=creal(out[ci]*conj(out[ci]));
				mb[1][ci]+=creal(out[ci]*conj(out[ci]));
			}	
		}
// First FFT (complex to complex)


// Second FFT (real to complex) à partir mb[1][ci] de  résultat dans mb[2][ci]
		for(ci=0;ci<NFFT;ci++) // copy SDR buffer to FFTW and no window (Vectorized)
				//ini[ci]=mb[1][ci];
				ini[ci]=sqrt(mb[1][ci]);// spectral power

		fftw_execute(pffti);

		for( ci=0;ci<NFFT;ci++) // vectorized
				mb[2][ci]=fabs(creal(outi[ci]*conj(outi[ci]))*creal(outi[ci]*conj(outi[ci])));
				//{
				//mb[2][ci]=creal(outi[ci]*conj(outi[ci]));
				//mb[3][ci]+=creal(outi[ci]*conj(outi[ci]));
				//}
// Second FFT (real to complex) à partir mb[1][ci] de  résultat dans mb[2][ci]

		for(ci=0;ci<NFFT;ci++) // freq matrix for pwl scaling
			fmat[cf][ci]=(frq+ci*RFstep-(frq_Samp/2))*1.0e-6; // MHz
		//for(ci=0;ci<NFFTd2;ci++) // Include FFT pos/neg swap
		//{
			//zrms[cf][ci+NFFTd2]=10*log10(mb[2][ci]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
			//zrms[cf][ci]=10*log10(mb[2][ci+NFFTd2]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
			//zrms[cf][ci+NFFTd2]=10*log10(mb[1][ci]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
			//zrms[cf][ci]=10*log10(mb[1][ci+NFFTd2]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
			//zmin[cf][ci+NFFTd2]=10*log10(mb[2][ci]+1e-20)+RNFFTdB-gaindB; // /NFFT
			//zmin[cf][ci]=10*log10(mb[2][ci+NFFTd2]+1e-20)+RNFFTdB-gaindB; // /NFFT
			//zmax[cf][ci+NFFTd2]=10*log10(mb[3][ci]+1e-20)+RNFFTdB-gaindB; // /NFFT
			//zmax[cf][ci]=10*log10(mb[3][ci+NFFTd2]+1e-20)+RNFFTdB-gaindB; // /NFFT
		//}
		for(ci=0;ci<NFFT;ci++) // Include FFT pos/neg swap
		{
			zrms[cf][ci]=10*log10(mb[2][ci]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT
//			zrms[cf][ci]=10*log10(mb[1][ci]+1e-20)+RNFFTdB+RNRptdB-gaindB; // /NFFT	
//			printf("%.0f\n",zrms[cf][ci]);
		}
		
// AFC calculation (end)  **************************************************************************************

		bord=500;

// calcul de la moyenne
// moyenne glissante sur 2*bord+1 points
// sur tous les points sauf les "bord" premiers et "bord" derniers

		for (ci=0;ci<NFFT;ci++)
			zmin[cf][ci]=0;
		for (ci=0;ci<bord;ci++)
			zmin[cf][ci]=zrms[cf][ci];
		for (ci=NFFT-bord;ci<NFFT;ci++)
			zmin[cf][ci]=zrms[cf][ci];
		for (ci=bord;ci<NFFT-bord;ci++)
		{
			for (ck=0;ck<(2*bord)+1;ck++)
				zmin[cf][ci]=zmin[cf][ci]+((zrms[cf][ci+ck-bord])/(2*bord+1));
		}
// calcul de la moyenne

// calcul ecart-type
		ecart=0;qecart=1.2;
		for(ci=bord;ci<NFFT-bord;ci++) 
			ecart=ecart+(zrms[cf][ci]-zmin[cf][ci])*(zrms[cf][ci]-zmin[cf][ci]);
		ecart=sqrt(ecart/NFFT);
// calcul ecart-type

// autocorrelation peak calculation (start)
// si valeur supérieur à la moyenne + X x ecart-type: c'est un pic sinon ce n'est aps un pic => 0

		for(ci=0;ci<((NFFT/2));ci++) 
		{
			if (zrms[cf][ci] > (zmin[cf][ci]+qecart*ecart))
				{
				zmax [cf][ci] = fabs(zrms[cf][ci]- zmin[cf][ci]);
				}		
			if (zrms[cf][ci] < (zmin[cf][ci]-qecart*ecart))
				{
				zmax [cf][ci] = fabs(zrms[cf][ci]- zmin[cf][ci]);
				zmax [cf][ci] = 0;
				}
			if ((zrms[cf][ci] >= (zmin [cf][ci]-qecart*ecart)) && (zrms[cf][ci] <= (zmin [cf][ci]+qecart*ecart)))
				zmax [cf][ci] = 0;
			if (Graph == 0)
				zrms[cf][ci]= zmax[cf][ci];
			if (Graph == 2)
				zrms[cf][ci]= zmin[cf][ci];
			if (Graph == 3)
				zrms[cf][ci]= zmax[cf][ci];
			
		}
// autocorrelation peak calculation (end)

// à enlever si mode Graphique ou tableau


// DMR identification (1) (start) ****************************************************************************

		if (Graph == 0)
			{
			//printf("Analysis : %i\n",cf+1);
			{
			struct tm *myTime=gmtime(&swpTime[0]); // utc time
		//	struct tm *myTime=localtime(&swpTime[0]);
			char dateStr[100];
			char time24Str[100]; // yr-mnth-dy,24h time
			sprintf(dateStr,"%i-%i-%i",myTime->tm_year+1900,myTime->tm_mon+1,myTime->tm_mday);
			sprintf(time24Str,"%i:%i:%i",myTime->tm_hour,myTime->tm_min,myTime->tm_sec);
			fprintf(stdout,"%s, %s, %.6e Hz: analysis n°%i\n",dateStr,time24Str,frq,cf+1);
			}
			compteur_DMR=0; seuil_DMR=6;
			if (zmin[cf][3000]-zmin[cf][4500] > seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 3000:  %.0f\n",fabs(zmin[cf][3000]-zmin[cf][4500]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+abs(zmin[cf][3000]-zmin[cf][4500])
				}
			if (zmin[cf][4500]-zmin[cf][6000] < -seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 4500:  %.0f\n",fabs(zmin[cf][4500]-zmin[cf][6000]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+fabs(zmin[cf][4500]-zmin[cf][6000])
				}
			if (zmin[cf][6000]-zmin[cf][7500] > seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 6000:  %.0f\n",fabs(zmin[cf][6000]-zmin[cf][7500]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+fabs(zmin[cf][6000]-zmin[cf][7500])
				}
			if (zmin[cf][7500]-zmin[cf][9000] < -seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 7500:  %.0f\n",fabs(zmin[cf][7500]-zmin[cf][9000]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+fabs(zmin[cf][7500]-zmin[cf][9000])
				}
			if (zmin[cf][9000]-zmin[cf][10500] > seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 9000:  %.0f\n",fabs(zmin[cf][9000]-zmin[cf][10500]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+fabs(zmin[cf][9000]-zmin[cf][10500])
				}
			if (zmin[cf][10500]-zmin[cf][12000] < -seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 10500:  %.0f\n",fabs(zmin[cf][10500]-zmin[cf][12000]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+fabs(zmin[cf][10500]-zmin[cf][12000])
				}
			if (zmin[cf][12000]-zmin[cf][13500] > seuil_DMR)
				{
				if (Verbose == 2)
					printf("DMR point 12000:  %.0f\n",fabs(zmin[cf][12000]-zmin[cf][13500]));
				compteur_DMR=compteur_DMR+1;
				//DMR_C=DMR_C+1;
				//DMR_L=DMR_L+fabs(zmin[cf][12000]-zmin[cf][13500])
				}
			//printf("Compteur DMR:    %i\n",compteur_DMR);
			if (compteur_DMR > 3)
				{
				if (Verbose == 2)
					printf("DMR recognized\n");
				DMR_C=DMR_C+1;
				DMR_L=0;
				}
			}

// DMR identification (1) (start) *************************************************************************************


// other identification then DMR (2) (start) **************************************************************************

		if (compteur_DMR < 4)
		{
		if (Graph == 0)
		{

//TETRAPOL identification - 8 tests (start) ****************************************************************************

		if (zrms[cf][1000]>0.0)
			{
			if (Verbose == 2)
				printf("TETRAPOL point 1000:  %.0f\n",zrms[cf][1000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][1000];
			}
		if (zrms[cf][2000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 2000:  %.0f\n",zrms[cf][2000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][2000];
			}	
		if (zrms[cf][10000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 10000:  %.0f\n",zrms[cf][10000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][10000];
			}
		if (zrms[cf][11000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 11000:  %.0f\n",zrms[cf][11000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][11000];
			}
		if (zrms[cf][4000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 4000:  %.0f\n",zrms[cf][4000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][4000];
			}
		if (zrms[cf][5000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 5000:  %.0f\n",zrms[cf][5000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][5000];
			}		
		if (zrms[cf][8000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 8000:  %.0f\n",zrms[cf][8000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][8000];
			}
		if (zrms[cf][7000]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRAPOL point 7000:  %.0f\n",zrms[cf][7000]);
			TETRAPOL_C= TETRAPOL_C+1;
			TETRAPOL_L= TETRAPOL_L+zrms[cf][7000];
			}

//TETRAPOL identification - 8 tests (start) ****************************************************************************

//MPT1237 identification - 2 tests (start) *****************************************************************************

		if (zrms[cf][5333]>0.0) 
			{
			if (Verbose == 2)
				printf("STANAG 4285  or MPT1237 point 5333:  %.0f\n",zrms[cf][5333]);
			MPT1237_C= MPT1237_C+1;
			MPT1237_L= MPT1237_L+zrms[cf][5333];
			}
		if (zrms[cf][10667]>0.0) 
			{
			if (Verbose == 2)
				printf("STANAG 4285  or MPT1237 point 10667:  %.0f\n",zrms[cf][10667]);
			MPT1237_C= MPT1237_C+1;
			MPT1237_L= MPT1237_L+zrms[cf][10667];
			}

//MPT1237 identification - 2 tests (end) ******************************************************************************

//TETRA identification - 7 tests (start) ******************************************************************************

		if (zrms[cf][2833]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 2833:  %.0f\n",zrms[cf][2833]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][2833];
			};
		if (zrms[cf][5667]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 5667:  %.0f\n",zrms[cf][5667]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][5667];
			};
		if (zrms[cf][8500]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 8500:  %.0f\n",zrms[cf][8500]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][8500];
			};
		if (zrms[cf][11333]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 11333:  %.0f\n",zrms[cf][11333]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][11333];
			};

		if (zrms[cf][6375]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 6375:  %.0f\n",zrms[cf][6375]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][6375];
			};
		if (zrms[cf][7083]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 7083:  %.0f\n",zrms[cf][7083]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][7083];
			};
		if (zrms[cf][7792]>0.0) 
			{
			if (Verbose == 2)
				printf("TETRA point 7792:  %.0f\n",zrms[cf][7792]);
			TETRA_C= TETRA_C+1;
			TETRA_L= TETRA_L+zrms[cf][7792];
			};

//TETRA identification - 7 tests (end) ******************************************************************************
		
//DMR identification - 4 tests (start) ******************************************************************************

/*
		//if (zrms[cf][12000]>0.0) 
			{
			printf("DMR point 12000:  %.0f\n",zmin[cf][12000]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][12000];
			};
		//if (zrms[cf][3000]>0.0) 
			{
			printf("DMR point 3000:  %.0f\n",zmin[cf][3000]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][3000];
			};		
		//if (zrms[cf][9000]>0.0) 
			{
			printf("DMR point 9000:  %.0f\n",zmin[cf][9000]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][9000];
			};
		//if (zrms[cf][6000]>0.0) 
			{
			printf("DMR point 6000:  %.0f\n",zmin[cf][6000]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][6000];
			};

		//if (zrms[cf][4500]>0.0) 
			{
			printf("DMR point 4500:  %.0f\n",zmin[cf][4500]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][12000];
			};
		//if (zrms[cf][7500]>0.0) 
			{
			printf("DMR point 7500:  %.0f\n",zmin[cf][7500]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][3000];
			};		
		//if (zrms[cf][10500]>0.0) 
			{
			printf("DMR point 10500:  %.0f\n",zmin[cf][10500]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][9000];
			};
		//if (zrms[cf][13500]>0.0) 
			{
			printf("DMR point 13500:  %.0f\n",zmin[cf][13500]);
			DMR_C= DMR_C+1;
			DMR_L= DMR_L+zrms[cf][6000];
			};
*/

//DMR identification - 4 tests (end) ******************************************************************************

//POCSAG identification - 4 tests (start) *************************************************************************

		if (zrms[cf][1333]>0.0) 
			{
			if (Verbose == 2)
				printf("POCSAG point 1333:  %.0f\n",zrms[cf][1333]);
			POCSAG_C= POCSAG_C+1;
			POCSAG_L= POCSAG_L+zrms[cf][1333];
			};
		if (zrms[cf][2667]>0.0) 
			{
			if (Verbose == 2)
				printf("POCSAG point 2667:  %.0f\n",zrms[cf][2667]);
			POCSAG_C= POCSAG_C+1;
			POCSAG_L= POCSAG_L+zrms[cf][2667];
			};
		if (zrms[cf][6667]>0.0) 
			{
			if (Verbose == 2)
				printf("POCSAG point 6667:  %.0f\n",zrms[cf][6667]);
			POCSAG_C= POCSAG_C+1;
			POCSAG_L= POCSAG_L+zrms[cf][6667];
			};
		if (zrms[cf][5333]>0.0) 
			{
			if (Verbose == 2)
				printf("POCSAG point 5333:  %.0f\n",zrms[cf][5333]);
			POCSAG_C= POCSAG_C+1;
			POCSAG_L= POCSAG_L+zrms[cf][5333];
			};

//POCSAG identification - 4 tests (end) *************************************************************************

//STANDARDC identification - 2 tests (start) ********************************************************************

		if (zrms[cf][13500]>0.0) 
			{
			if (Verbose == 2)
				printf("StandardC point 13500:  %.0f\n",zrms[cf][13500]);
			STANDARDC_C= STANDARDC_C+1;
			STANDARDC_L= STANDARDC_L+zrms[cf][13500];
			};
		if (zrms[cf][6750]>0.0) 
			{
			if (Verbose == 2)
				printf("StandardC point 6750:  %.0f\n",zrms[cf][6750]);
			STANDARDC_C= STANDARDC_C+1;
			STANDARDC_L= STANDARDC_L+zrms[cf][6750];
			}

//STANDARDC identification - 2 tests (end) *********************************************************************

//GSM identification - 4 tests (start) ***************************************************************************

		if (zrms[cf][4615]>0.0) 
			{
			if (Verbose == 2)
				printf("GSM point 4615:  %.0f\n",zrms[cf][4615]);
			GSM_C= GSM_C+1;
			GSM_L= GSM_L+zrms[cf][4615];
			}
		if (zrms[cf][3692]>0.0) 
			{
			if (Verbose == 2)
				printf("GSM point 3692:  %.0f\n",zrms[cf][3692]);
			GSM_C= GSM_C+1;
			GSM_L= GSM_L+zrms[cf][3692];
			}
		if (zrms[cf][2307]>0.0) 
			{
			if (Verbose == 2)
				printf("GSM point 2307:  %.0f\n",zrms[cf][2307]);
			GSM_C= GSM_C+1;
			GSM_L= GSM_L+zrms[cf][2307];
			}
		if (zrms[cf][1846]>0.0) 
			{
			if (Verbose == 2)
				printf("GSM point 1846:  %.0f\n",zrms[cf][1846]);
			GSM_C= GSM_C+1;
			GSM_L= GSM_L+zrms[cf][1846];
			}

//GSM identification - 4 tests (end) ***************************************************************************

//DAB identification - 3 tests (start) ***************************************************************************

		if (zrms[cf][4800]>0.0) 
			{
			if (Verbose == 2)
				printf("DAB point 4800:  %.0f\n",zrms[cf][4800]);
			DAB_C= DAB_C+1;
			DAB_L= DAB_L+zrms[cf][4800];
			}
		if (zrms[cf][9600]>0.0) 
			{
			if (Verbose == 2)
				printf("DAB point 9600:  %.0f\n",zrms[cf][9600]);
			DAB_C= DAB_C+1;
			DAB_L= DAB_L+zrms[cf][9600];
			}
		if (zrms[cf][14400]>0.0) 
			{
			if (Verbose == 2)
				printf("DAB point 14400:  %.0f\n",zrms[cf][14400]);
			DAB_C= DAB_C+1;
			DAB_L= DAB_L+zrms[cf][14400];
			}

//DAB identification - 3 tests (end) ***************************************************************************

//LTE identification - 8 tests (start) ***************************************************************************

		if (zrms[cf][1500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 1500:  %.0f\n",zrms[cf][1500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][1500];
			}
		if (zrms[cf][2500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 2500:  %.0f\n",zrms[cf][2500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][2500];
			}
		if (zrms[cf][3500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 3500:  %.0f\n",zrms[cf][3500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][3500];
			}
		if (zrms[cf][4500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 4500:  %.0f\n",zrms[cf][4500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][4500];
			}
		if (zrms[cf][5500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 5500:  %.0f\n",zrms[cf][5500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][5500];
			}
		if (zrms[cf][6500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 6500:  %.0f\n",zrms[cf][6500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][6500];
			}
		if (zrms[cf][7500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 7500:  %.0f\n",zrms[cf][7500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][7500];
			}
		if (zrms[cf][8500]>0.0)
			{
			if (Verbose == 2)
				printf("LTE point 8500:  %.0f\n",zrms[cf][8500]);
			LTE_C= LTE_C+1;
			LTE_L= LTE_L+zrms[cf][8500];
			}

//LTE identification - 8 tests (end) *****************************************************************************

		}
		}

// other identification then DMR (2) (stop) *********************************************************************

		if (Graph == 0)
			printf("******************************************************************************\n");

	}

// AFC Analysis (B) (end) ****************************************************************************************

	
//	for( ci=0;ci<4;ci++)
//		free(mb[ci]);
//	free(buf);
	
	fCnt=fCnt+1; //historic reason


//Modulation identification -final (3) (start) ********************************************************************	

	if (Graph == 0)
	{
		if (Verbose != 0)
			{
			printf("TETRAPOL       : %i %.0f\n", (100*TETRAPOL_C/((fCnt-1)*8)), TETRAPOL_L/TETRAPOL_C);
			printf("TETRA          : %i %.0f\n", (100*TETRA_C/((fCnt-1)*7)), TETRA_L/TETRA_C);
			printf("DMR            : %i %.0f\n", (100*DMR_C/((fCnt-1)*1)), DMR_L/DMR_C);
			printf("MPT1237        : %i %.0f\n", (100*MPT1237_C/((fCnt-1)*2)), MPT1237_L/MPT1237_C);
			printf("POCSAG         : %i %.0f\n", (100*POCSAG_C/((fCnt-1)*4)), POCSAG_L/POCSAG_C);
			printf("STANDARDC      : %i %.0f\n", (100*STANDARDC_C/((fCnt-1)*2)), STANDARDC_L/STANDARDC_C);
			printf("GSM            : %i %.0f\n", (100*GSM_C/((fCnt-1)*4)), GSM_L/GSM_C);
			printf("DAB            : %i %.0f\n", (100*DAB_C/((fCnt-1)*3)), DAB_L/DAB_C);
			printf("LTE            : %i %.0f\n", (100*LTE_C/((fCnt-1)*8)), LTE_L/LTE_C);
			printf("******************************************************************************\n");
			}
	

	//max= 100*TETRAPOL_C/((fCnt-1)*8); winner = 1;
	max= 0; winner = 0;
	if (((100*TETRAPOL_C/((fCnt-1)*8)) > max) && ((100*LTE_C/((fCnt-1)*8)) < 49))
		{
		max= 100*TETRAPOL_C/((fCnt-1)*8); winner = 1;
		}
// new V4 (start)
	if (((100*TETRAPOL_C/((fCnt-1)*8)) > max) && ((100*LTE_C/((fCnt-1)*8)) >= 49))
		{
		if ((TETRAPOL_L/TETRAPOL_C) > (1.4*(LTE_L/LTE_C)))
			{
			max= 100*TETRAPOL_C/((fCnt-1)*8); winner = 1;
			}
		}
// new V4 (end)
	if ((100*TETRA_C/((fCnt-1)*7)) > max)
		{
		max= 100*TETRA_C/((fCnt-1)*7); winner = 2;
		}
	if ((100*DMR_C/((fCnt-1)*1)) > max)
		{
		max= 100*DMR_C/((fCnt-1)*1); winner = 3;
		}
	if (((100*MPT1237_C/((fCnt-1)*2)) > max) && ((100*POCSAG_C/((fCnt-1)*4)) < 50))
		{
		max= 100*MPT1237_C/((fCnt-1)*2); winner = 4;
		}
	if ((100*POCSAG_C/((fCnt-1)*4)) > max)
		{
		max= 100*POCSAG_C/((fCnt-1)*4); winner = 5;
		}
	if (((100*STANDARDC_C/((fCnt-1)*2)) > max) && ((100*GSM_C/((fCnt-1)*4)) < 50))
		{
		max= 100*STANDARDC_C/((fCnt-1)*2); winner = 6;
		}
	if ((100*GSM_C/((fCnt-1)*4)) > max)
		{
		max= 100*GSM_C/((fCnt-1)*4); winner = 7;
		}
	if ((100*DAB_C/((fCnt-1)*3)) > max)
		{
//		if ((DAB_L/DAB_C) > 15)
			max= 100*DAB_C/((fCnt-1)*3); winner = 8;
		}
	if ((100*LTE_C/((fCnt-1)*8)) > max) 
		{
		max= 100*LTE_C/((fCnt-1)*8); winner = 9;
		}

	if (max < 51)
		winner = 0;

	//printf("%.0f %.0f\n",max,winner);
	{
	struct tm *myTime=gmtime(&swpTime[0]); // utc time
//	struct tm *myTime=localtime(&swpTime[0]);
	char dateStr[100];
	char time24Str[100]; // yr-mnth-dy,24h time
	sprintf(dateStr,"%i-%i-%i",myTime->tm_year+1900,myTime->tm_mon+1,myTime->tm_mday);
	sprintf(time24Str,"%i:%i:%i",myTime->tm_hour,myTime->tm_min,myTime->tm_sec);
	fprintf(stdout,"%s, %s, %.6e Hz : ",dateStr,time24Str,frq);
	}
	if (winner == 0)
		printf("\033[31;01m No signal recognized\033[00m\n");
	if (winner == 1)
		printf("\033[32;01m TETRAPOL signal recognized\033[00m\n");
	if (winner == 2)
		printf("\033[32;01m TETRA signal recognized\033[00m\n");
	if (winner == 3)
		printf("\033[32;01m DMR signal recognized\033[00m\n");
	if (winner == 4)
		{
		if (frq > 30000)
			printf("\033[32;01m MPT1237 signal recognized\033[00m\n");
		if (frq <= 30000)
			printf("\033[32;01m STANAG4285 signal recognized\033[00m\n");
		}
	if (winner == 5)
		printf("\033[32;01m POCSAG signal recognized\033[00m\n");
	if (winner == 6)
		printf("\033[32;01m STANDARDC signal recognized\033[00m\n");
	if (winner == 7)
		printf("\033[32;01m GSM signal recognized\033[00m\n");
	if (winner == 8)
		printf("\033[32;01m DAB signal recognized\033[00m\n");
	if (winner == 9)
		printf("\033[32;01m LTE signal recognized\033[00m\n");
	printf("******************************************************************************\n");
	}

//Modulation identification -final (3) (end) *******************************************************************
	
	fCnt=fCnt-1; // historic reason
	
    	TETRAPOL_C=0, TETRAPOL_L=0;
	TETRA_C=0, TETRA_L=0;
	DMR_C=0, DMR_L=0;
	MPT1237_C=0, MPT1237_L=0;
	POCSAG_C=0, POCSAG_L=0;
	STANDARDC_C=0, STANDARDC_L=0;
	DVBT_C=0, DVBT_L=0;
	GSM_C=0, GSM_L=0;
	DAB_C=0, DAB_L=0;
	LTE_C=0, LTE_L=0;
	compteur_DMR=0; seuil_DMR=6;
	
	if (NScan != 1)
		winner = 99;
	} while (winner == 0);

// loop do-while (A) (end)**************************************************************************************

	for( ci=0;ci<4;ci++)
		free(mb[ci]);
	free(buf);

}

void GnuPlotDispAndSave( char *fName,double **zz )
{ // zz[NY][NX], fName != fNameStem, as we have several sets of files to plot/print
	unsigned int ci,cj;
	char fname[100];
	FILE *ftdv, *fcsv;
	float fftsf=(frq_Samp*1.0e-6)/NFFT;
	float frq=frq_min;
	FILE *GpPipe;
#ifdef USE_GNUPLOT
//-	GpPipe=popen("gnuplot","w"); // init GnuPlot Pipe
//-	fprintf(GpPipe,"set term gif\n");
//-	fprintf(GpPipe,"set output \"%s/%s.gif\"\n",fNameStem,fName);
//-	fprintf(GpPipe,"set grid\n");
//-	fprintf(GpPipe,"set surface\n"); // not required?
//	fprintf(GpPipe,"set view 10,340\n");
//	fprintf(GpPipe,"set xrange [0:10]\n");
//	fprintf(GpPipe,"set yrange [0:10]\n");
//	fprintf(GpPipe,"set zrange [-90:0]\n");
//-	fprintf(GpPipe,"set pm3d\n"); // add map for heat map
//-	fprintf(GpPipe,"set xlabel \"FFT Band MHz\"\n");
//-	fprintf(GpPipe,"set ylabel \"FFT Center MHz\"\n");
//-	fprintf(GpPipe,"set xtics 2\n");
//-	fprintf(GpPipe,"set ytics %f\n",frq_Tic);
//	fprintf(GpPipe,"set palette defined (-1 \"blue\", 0 \"green\", 1 \"yellow\", 2 \"red\")\n");
//	fprintf(GpPipe,"set palette defined (-3 \"dark-green\", -2 \"green\", -1 \"cyan\", 0 \"blue\", 1 \"purple\", 2 \"red\", 3 \"orange\", 4 \"white\")\n");
//-	fprintf(GpPipe,"set palette defined (-5 \"black\",-4 \"dark-green\",-3 \"forest-green\", -2 \"green\", -1 \"sea-green\", 0 \"cyan\", 1 \"blue\", 2 \"purple\", 3 \"pink\",4 \"red\", 5 \"orange\", 6 \"yellow\", 7 \"white\")\n");
//	fprintf(GpPipe,"set palette defined ( 0 0 0 0, 1 1 1 1 )\n"); // grey scale
//	fprintf(GpPipe,"set palette defined ( 0 0 0 0, 1 0 1 1 )\n"); // cyan scale
//-	fprintf(GpPipe,"set hidden3d\n");
//	fprintf(GpPipe,"set contour base\n"); // base surface both
//	fprintf(GpPipe,"splot '-' using 1:2:3 with pm3d title 'ch 1'\n");
//	for(ci=0;ci<fCnt;++ci)
//	{ // 1-D x y z format
//		for(cj=0;cj<NFFT;++cj)
//			fprintf(GpPipe,"%f %f %f\n",(1.0*cj-(NFFT/2))*fftsf,frq*1.0e-6,zz[ci][cj]);
//		fprintf(GpPipe,"\n");
//		frq+=frq_step;
//	}
//	fprintf(GpPipe,"e\n");
//	fflush(GpPipe);
//	pclose(GpPipe); // kill gnuplot process!
#endif
//	sprintf( fname, "%s/%s.xls",fNameStem,fName);
//
//	ftdv=fopen(fname,"w"); // standard xls tab delimited variable "\t" format
//	fprintf(ftdv,"\t\t");
//	for(cj=0;cj<NFFT;++cj)
//		fprintf(ftdv,"\t%i",cj);
//	fprintf(ftdv,"\n");
//	for(ci=0;ci<fCnt;++ci)
//	{
//		fprintf(ftdv,"%.3f\t%.3f\t%i\t",(frq_min+frq_step*(ci-0.5))/1e6,(frq_min+frq_step*(ci+0.5))/1e6,ci);
//		for(cj=0;cj<NFFT;++cj)
//			fprintf(ftdv,"%.2f\t",zz[ci][cj]);
//		fprintf(ftdv,"\n");
//	}
//	fclose(ftdv);

//	sprintf( fname, "%s/%s.csv",fNameStem,fName);
	if( strcmp( fNameStem,"output")!=0 )
		fcsv=fopen(fname,"w"); // "soapy power" style csv ", " format
	else
		fcsv=stdout;
//	fprintf(fcsv,"%s\n",fname);

// changement fCnt
//	fCnt=fCnt-1;

	for(ci=0;ci<fCnt;++ci)
	{
		struct tm *myTime=gmtime(&swpTime[ci]); // utc time
//		struct tm *myTime=localtime(&swpTime[ci]);
		char dateStr[100];
		char time24Str[100]; // yr-mnth-dy,24h time
		sprintf(dateStr,"%i-%i-%i",myTime->tm_year+1900,myTime->tm_mon+1,myTime->tm_mday);
		sprintf(time24Str,"%i:%i:%i",myTime->tm_hour,myTime->tm_min,myTime->tm_sec);
		if (Graph != 0)
			fprintf(fcsv,"%s, %s, ",dateStr,time24Str);

//ligne à décommenter en mode tableau ou Graphique
		if (Graph != 0)
			fprintf(fcsv,"%.0f, %.0f, %.1f, %i",(frq_min+frq_step*(ci-0.5)),(frq_min+frq_step*(ci+0.5)),frq_Samp/NFFT,NFFTcro); // frq_min, frq_max, fftbin_frq_width, buffer_length_samples
//ligne à décommenter en mode tableau ou Graphique

		//fprintf(fcsv,"%.3f, %.3f, %.3f, %i",(0.0,(10e6*1000.0*NFFT/2*frq_Samp),(10e6*1000.0/frq_Samp),NFFTcro)); // frq_min, frq_max, fftbin_frq_width, buffer_length_samples
		//fprintf(fcsv,"%.0f, %.0f, %.1f, %i",(655.0,1310.0,0.02,NFFTcro)); // frq_min, frq_max, fftbin_frq_width, buffer_length_samples
		for(cj=(NFFT-NFFTcro)/2;cj<(NFFT/2)-6;++cj)
		{
			//if ((cj!=0) && (cj!=NFFT))
				//{
					//if ((zz[ci][cj]-zz[ci][cj-1] > 4.0) && (zz[ci][cj]-zz[ci][cj+1]) > -4.0)
						//{
							//printf("spuriousA%.2f\n",zz[ci][cj]-zz[ci][cj-1]);
							//printf("spuriousB%.2f\n",zz[ci][cj]-zz[ci][cj+1]);
							//zz[ci][cj] = (zz[ci][cj-2]+zz[ci][cj+2])/2;
						//}
					//if ((zz[ci][cj]-zz[ci][cj-2] > 6.0) && (zz[ci][cj]-zz[ci][cj+2]) > -6.0)
						//{
							//printf("spuriousA%.2f\n",zz[ci][cj]-zz[ci][cj-1]);
							//printf("spuriousB%.2f\n",zz[ci][cj]-zz[ci][cj+1]);
							//zz[ci][cj] = (zz[ci][cj-2]+zz[ci][cj+2])/2;
							//zz[ci][cj-1] = (zz[ci][cj-2]+zz[ci][cj+2])/2;
							//zz[ci][cj+1] = (zz[ci][cj-2]+zz[ci][cj+2])/2;
						//}
					//if ((zz[ci][cj]-zz[ci][cj-4] > 4.0) && (zz[ci][cj]-zz[ci][cj+4]) > -4.0)
						//{
							//printf("spuriousA%.2f\n",zz[ci][cj]-zz[ci][cj-1]);
							//printf("spuriousB%.2f\n",zz[ci][cj]-zz[ci][cj+1]);
							//zz[ci][cj] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
							//zz[ci][cj-1] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
							//zz[ci][cj+1] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
							//zz[ci][cj-2] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
							//zz[ci][cj+2] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
							//zz[ci][cj-3] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
							//zz[ci][cj+3] = (zz[ci][cj-4]+zz[ci][cj+4])/2;
						//}
				//}

//ligne à décommenter en mode tableau ou Graphique
			if (Graph != 0)
				fprintf(fcsv,", %.2f",zz[ci][cj]);
//ligne à décommenter en mode tableau ou Graphique
		}
		if (Graph != 0)
			fprintf(fcsv,"\n");
	}
	if( fcsv!=stdout )
		fclose(fcsv);
}

void DecCmdLine( int argc, char *argv[] )
{ // support a subset of Soapy Power for compatibility
	int ci=0;
	char units1,units2,sep;
	int scanLen;
	char fname[100];
	signed char quitOnHelp=0;
	FILE *fp;
	while( ci<argc )
	{
		if( strcmp(argv[ci],"-A" )==0 ) // antenna/LNA name
		{
			ci++;
			if( strcmp(argv[ci],"LNAH" )==0 )
				LNAnum=1;
			if( strcmp(argv[ci],"LNAL" )==0 )
				LNAnum=2;
			if( strcmp(argv[ci],"LNAW" )==0 )
				LNAnum=3;
		}
//		if( strcmp(argv[ci],"-b" )==0 ) // number of FFT bins
//		{
//			ci++;
//			NFFT=atoi(argv[ci]);
//		}
		if( strcmp(argv[ci],"-verb" )==0 ) // Verbose = 0 (no) - Verbose = 1 (low) - Verbose = 2 (moderate)
		{
			ci++;
			Verbose=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-C" )==0 ) // channel number - always 0 for LimeSDRmini
		{
			ci++;
			Ch=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-n" )==0 ) // repeats
		{
			ci++;
			NRpt=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-ns" )==0 ) // scan
		{
			ci++;
			NScan=atoi(argv[ci]);
		}
                if( strcmp(argv[ci],"-sg" )==0 ) // mesure scan speed
		{
			ci++;
			StopGo=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-d" )==0 ) // scan
		{
			ci++;
			Dev=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-w" )==0 ) // RF filter width Hz
		{
			ci++;
			scanLen=sscanf(argv[ci],"%f%c",&frq_LPF,&units1);
			if( scanLen>1 ) // detect if units used
			{
				if( units1=='k' )
					frq_LPF*=1e3;
				if( units1=='M' )
					frq_LPF*=1e6;
				if( units1=='G' )
					frq_LPF*=1e9;
			}
		}
		if( strcmp(argv[ci],"-OSR" )==0 ) // Oversample rate
		{
			ci++;
			OSR=atoi(argv[ci]);
		}
		if( strcmp(argv[ci],"-st" )==0 ) // Scan step
		{
			ci++;
			scanLen=sscanf(argv[ci],"%f%c",&scan_step,&units1);
			if( scanLen>1 ) // detect ifunits used
			{
				if( units1=='k' )
					scan_step*=1e3;
				if( units1=='M' )
					scan_step*=1e6;
				if( units1=='G' )
					scan_step*=1e9;
			}
		}
//		if( strcmp(argv[ci],"-r" )==0 ) // Sample Rate at USB S/s
//		{
//			ci++;
//			scanLen=sscanf(argv[ci],"%f%c",&frq_Samp,&units1);
//			if( scanLen>1 ) // detect ifunits used
//			{
//				if( units1=='k' )
//					frq_Samp*=1e3;
//				if( units1=='M' )
//					frq_Samp*=1e6;
//				if( units1=='G' )
//					frq_Samp*=1e9;
//			}
//		}
		if( strcmp(argv[ci],"-g" )==0 ) // Gain dB
		{
			ci++;
			gaindB=atoi(argv[ci]);
		}
		if( (strcmp(argv[ci],"-O" )==0) || (strcmp(argv[ci],"-o" )==0) ) // Output Filename
		{
			ci++;
			strcpy(fNameStem,argv[ci]);
		}
		if( strcmp(argv[ci],"-f" )==0 ) // frequency range, FORTRAN style fmin:fmax
		{ // sopay power also supports a center frequency mode, we do not!
			ci++;
			units1=' ';
			units2=' ';
			scanLen=sscanf(argv[ci],"%f%c%c%f%c",&frq_min,&units1,&sep,&frq_max,&units2);
			if( scanLen>3 ) // detect if units used
			{
				if( units1=='k' )
					frq_min*=1.0e3;
				if( units1=='M' )
					frq_min*=1.0e6;
				if( units1=='G' )
					frq_min*=1.0e9;
				if( units2=='k' )
					frq_max*=1.0e3;
				if( units2=='M' )
					frq_max*=1.0e6;
				if( units2=='G' )
					frq_max*=1.0e9;
				//printf("%0;f\n",frq_min);  il y a 8 Hz de plus????
			}
		}
		if( (strcmp(argv[ci],"-v" )==0) || (strcmp(argv[ci],"--version" )==0) )
			printf("Version=%f\n",0.0); // Version number
		if( (strcmp(argv[ci],"-h" )==0) || (strcmp(argv[ci],"--help" )==0) )
		{
			printf("LimeScan [-h] [--help] [-f Hz:Hz] [-O FILEstub] [-o FILEstub] [--info] [-v] [--version] [-b BINS] [-w Hz] [-r S/s] [-n REPEATS] [-t SECONDS] [-A ANTENNA] [-C CHANNEL]\n");
			quitOnHelp=1;
		}
		if( strcmp(argv[ci],"--info" )==0 ) // Soapy device info - not relevant
			printf("No information available, please consult your system administrator\n");
		if( strcmp(argv[ci],"-gps" )==0 )
		{
			ci++;
			strcpy(gpsPath,argv[ci]);
			doGPS=1;
		}
		if( strcmp(argv[ci],"-Tst" )==0 )
		{
			ci++;
			tstLvl=atoi(argv[ci]);
			doTst=1;
			printf("-Tst %i\n",tstLvl);
		}
		if( strcmp(argv[ci],"-pwla" )==0 )
			doPwlAnt=1;
		if( strcmp(argv[ci],"-pwlw" )==0 )
			doPwlLNAW=1;
		if( strcmp(argv[ci],"-pwlh" )==0 )
			doPwlLNAH=1;
		if( strcmp(argv[ci],"-pwll" )==0 )
			doPwlLNAL=1;
		ci++;
	}
	mkdir( fNameStem, S_IRWXU ); // S_IROTH Note: default ./output - for graphics and .xls; .csv to stdout
	sprintf(fname,"%s/%s_cmd_params.txt",fNameStem,fNameStem);
	fp=fopen(fname,"w");
	for(ci=0;ci<argc;ci++) // make copy of command line with data.
		fprintf(fp,"%s ",argv[ci]);
	fprintf(fp,"\n\n");
	printf("File=%s\n",fNameStem);
	if( LNAnum==1 )
		fprintf(fp,"LNAH Ch=%i\n",Ch);
	if( LNAnum==2 )
		fprintf(fp,"LNAL Ch=%i\n",Ch);
	if( LNAnum==3 )
		fprintf(fp,"LNAW Ch=%i\n",Ch);
	fprintf(fp,"Freq Range %.3f:%.3f MHz\n", frq_min*1.0e-6,frq_max*1.0e-6 );
	if(((frq_max-frq_min)/1.0E6)>=10)
		frq_Tic=2.0;
	if(((frq_max-frq_min)/1.0E6)>=30)
		frq_Tic=10.0;
	if(((frq_max-frq_min)/1.0E6)>=100)
		frq_Tic=20.0;
	if(((frq_max-frq_min)/1.0E6)>=300)
		frq_Tic=100.0;
	if(((frq_max-frq_min)/1.0E6)>=1000)
		frq_Tic=200.0;
	if(((frq_max-frq_min)/1.0E6)>=2000)
		frq_Tic=500.0;
	fprintf(fp,"Rate=%.3f Ms/s OSR=%i\n",frq_Samp*1.0e-6,OSR);
	fprintf(fp,"%f\n",4*OSR*(frq_Samp*1.0e-6));
	if( (4*OSR*(frq_Samp*1.0e-6))>640 )
		fprintf(fp,"ERROR: Rate*OSR must be < 640MS/s\n");
	fprintf(fp,"Gain=%idB\n",gaindB);
	fprintf(fp,"LPF=%.3f MHz (RF)\n",frq_LPF*1.0e-6);
	fprintf(fp,"FFTW %i bins, with Hann Window, NRpt=%i\n",NFFT,NRpt);
	fclose(fp);

	if(quitOnHelp>0)
		exit(0);
}

void DisplayBinFile( char *fname )
{ // xdg-open can open txt,pdf,png,jpg,gif,mpg,avi. Preinstalled in Ubuntu 16.04
	char cmd[255];
	int err;
	sprintf(cmd,"xdg-open %s",fname);
	printf("SYSTEM(\"%s\")\n",cmd);
	err=system(cmd); // alternatives include: feh,eog,shotwell,tycat,tiv,gnome-open,pxl
	if( err!=0 )
		printf("system() %i %s\n",err,strerror(errno));
}

int ReadPwl( char fname[],double xlim[],double mPwl[],double cPwl[] )
{ // file format 'freq MHz', 'level dB'
	double ylim[256];
	short ci;
	unsigned char limCnt=0;
	FILE *fp=NULL;
	if( (fp=fopen(fname,"r"))==NULL )
	{
		printf("PWL file %s not found\n",fname);
		exit(1);
	}
	while( 1 )
	{
		if( fscanf(fp,"%lf %lf",&xlim[limCnt],&ylim[limCnt])<2 )
			break;
		limCnt++;
	}
	if( limCnt<2 )
	{
		printf("ERROR - too few limits\n");
		exit(1);
	}
	for(ci=0;ci<(limCnt-1);ci++)
	{
		mPwl[ci]=(ylim[ci+1]-ylim[ci])/(xlim[ci+1]-xlim[ci]);
		cPwl[ci]=ylim[ci]-xlim[ci]*mPwl[ci];
	}
	return(limCnt-1);
}

void ApplyPwl( double *vecy[],double xlim[],double mPwl[],double cPwl[],unsigned char nPwl )
{ // linear interpolation between pwl table
	int ci;
	int cj;
	unsigned char ck=0;
	for(cj=0;cj<fCnt;cj++)
		for(ci=0;ci<NFFT;ci++)
		{
			while( (xlim[ck+1]<fmat[cj][ci]) && (ck<nPwl) ) // check in interval range
				ck++;
			vecy[cj][ci]-=fmat[cj][ci]*mPwl[ck]+cPwl[ck];
		}
}

void ReadGPRMC( char *buffer )
{ // time UTC, status, lat, long, speed, track, date,
	float gpsTime[3]={0,0,0}; // 0 hr, 1 min, 2 sec
	int gpsDate[3]={0,0,0}; // 0 day, 1 mnth, 2 yr
	short gpsLatDeg=0; // deg
	float gpsLatMin=0; // min
	char gpsLatRef=' ';
	short gpsLngDeg=0; // deg
	float gpsLngMin=0; // min
	char gpsLngRef=' ';
	float gpsTrackT=0; // heading relative to true north
	float gpsSpeedK=0; // km/h

	char *buf=buffer;
	char word[20];
	char subWord[20];
	if( buf!=NULL )
	{ // hh mm ss.ss
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsTime[0]=atof(subWord);
		ReadSubWord(word,subWord,2,2);
		gpsTime[1]=atof(subWord);
		ReadSubWord(word,subWord,4,5);
		gpsTime[2]=atof(subWord);
	}
	if( buf!=NULL )
	{ // hh mm ss.ss
		buf=ReadTil( buf, word );
		if( word[0]=='V' )
			printf("WARNING\n");
	}
	if( buf!=NULL )
	{ // latitude
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsLatDeg=atoi(subWord);
		ReadSubWord(word,subWord,2,10);
		gpsLatMin=atof(subWord);
	}
	if( buf!=NULL )
	{ // N or S
		buf=ReadTil( buf, word );
		gpsLatRef=word[0];
	}
	if( buf!=NULL )
	{ // long
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsLngDeg=atoi(subWord);
		ReadSubWord(word,subWord,2,10);
		gpsLngMin=atof(subWord);
	}
	if( buf!=NULL )
	{ // E or W
		buf=ReadTil( buf, word );
		gpsLngRef=word[0];
	}
	if( buf!=NULL )
	{ // speed
		buf=ReadTil( buf, word );
		gpsSpeedK=atof(word); // knots
	}
	if( buf!=NULL )
	{ // track
		buf=ReadTil( buf, word );
		gpsTrackT=atof(word);
	}
	if( buf!=NULL )
	{ // data
		buf=ReadTil( buf, word );
		ReadSubWord(word,subWord,0,2);
		gpsDate[0]=atoi(subWord);
		ReadSubWord(word,subWord,2,2);
		gpsDate[1]=atoi(subWord);
		ReadSubWord(word,subWord,4,2);
		gpsDate[2]=atoi(subWord);
	}
	printf("GPRMC %.0f:%.0f:%.2f UTC ",gpsTime[0],gpsTime[1],gpsTime[2]);
	printf("%i/%i/%i UK ",gpsDate[0],gpsDate[1],gpsDate[2]);
	printf("%io %.3f'%c %io %.3f'%c ",gpsLatDeg,gpsLatMin,gpsLatRef, gpsLngDeg,gpsLngMin,gpsLngRef);
	printf("vel=%.1fkt Hdg=%.3foT\n",gpsSpeedK,gpsTrackT);
}

char* ReadTil( char buffer[], char word[] )
{ // read til comma, newline or /0, occasionally get ,, * is start of checksum
	unsigned char ci=0;
	char *buf=buffer;
	word[0]='\0';
	while( ((*buf)!=',') && ((*buf)!='\n') && ((*buf)!='\0') && ((*buf)!='*') )
		word[ci++]=(*buf++);
	word[ci]='\0';
	if( ((*buf)=='\n') || ((*buf)=='\0') || ((*buf)=='*') )
		buf=NULL;
	else
		buf++;
	return(buf);
}

void ReadSubWord(char word[],char subWord[],int pos,int len)
{ // read len chars, or until end of string
	int ci=0;
	while( (ci<len) && (word[pos+ci]!='\0') )
	{
		subWord[ci]=word[pos+ci];
		ci++;
	}
	subWord[ci]='\0';
}

void ReadNMEAtraffic( char *gpsPath )
{
	int ci=0;
	int cj=0;
	int ck=0;
	char buf[256];
	size_t size;
	char byte='\0';
	char header[10];
	int fd=0;
	if( (fd=open(gpsPath, O_RDWR | O_NOCTTY | O_NDELAY ))<0 )
	{
		printf("Error - cannot open Ublox GPS %s\n",gpsPath );
		exit(1);
	}
	fcntl(fd,F_SETFL,0);
	for(ci=0; ci<150;ci++)
	{
		ck=0;
		byte='\0';
		while(byte!='\n')
		{
			while( (size=read(fd,&byte,1))<1 ); // wait til ready
			buf[ck++]=byte;
		}
		ck--;
		byte='\0';
		buf[ck]='\0';
		if( ck>0 ) // ignore empty lines, seems to send \n\n
		{
			for(cj=0;cj<6;cj++)
				header[cj]=buf[cj];
			header[6]='\0';
			if(strcmp(header,"$GPRMC")==0)
				ReadGPRMC( buf+7 );
		}
	}
	close(fd);
}

