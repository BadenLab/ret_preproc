#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_ParameterTable()

// 1 // check for Parameter Table
if (waveexists($"OS_Parameters")==1)
	print "OS_Parameters table already exists - just opening that..."
else
	print "Generating new OS_Parameters table..."

// make a new table
make /o/n=100 OS_Parameters = NaN

// reads data-header
wave wParamsNum
wave /T wParamsStr

// Define Entries
variable entry_position = 0

//////////// from the wParamsNum wave // Andre 2016 04 14
Variable xPixelsInd,yPixelsInd,realPixDurInd,lineDur,sampRate,sampPeriod,zoomIndx, pcName,delay,year,month,day,recday// Andre 2016 04 14 & 04 26
string recDate
string setup
/// GENERAL ////////////////////////////////////////////////////////////////////////////////////////////////

setdimlabel 0,entry_position,LineDuration,OS_Parameters
xPixelsInd = FindDimLabel(wParamsNum,0,"User_dxPix" )// Andre 2016 04 14
yPixelsInd = FindDimLabel(wParamsNum,0,"User_dyPix" )// Andre 2016 04 14
realPixDurInd = FindDimLabel(wParamsNum,0,"RealPixDur" )// Andre 2016 04 14
lineDur = (wParamsNum[xPixelsInd] *  wParamsNum[realPixDurInd]) * 10^-6// Andre 2016 04 14
OS_Parameters[%LineDuration] = lineDur
entry_position+=1

setdimlabel 0,entry_position,'samp_period',OS_Parameters
sampPeriod = (lineDur* wParamsNum[yPixelsInd])// Andre 2016 04 14
OS_Parameters[%samp_period] = sampPeriod
entry_position+=1

setdimlabel 0,entry_position,samp_rate_Hz,OS_Parameters
sampRate = 1/sampPeriod// Andre 2016 04 14
OS_Parameters[%samp_rate_Hz] = sampRate
entry_position+=1

SetDimLabel 0,entry_position,Data_Channel,OS_Parameters
OS_Parameters[%Data_Channel] = 0 // Fluorescence Data in wDataChX - default 0
entry_position+=1

SetDimLabel 0,entry_position,Data_Channel2,OS_Parameters
OS_Parameters[%Data_Channel2] = 1 // Fluorescence Data in wDataChX - default 1 (for Ratiometric)
entry_position+=1

SetDimLabel 0,entry_position,Trigger_Channel,OS_Parameters
OS_Parameters[%Trigger_Channel] = 2 // Trigger Data in wDataChX - default 2 
entry_position+=1

SetDimLabel 0,entry_position,Display_Stuff,OS_Parameters
OS_Parameters[%Display_Stuff] = 1 // generate graphs? - 0/1 - default 1
entry_position+=1

SetDimLabel 0,entry_position,LightArtifact_cut,OS_Parameters
OS_Parameters[%LightArtifact_cut] = 5 // nPixels cut in X to remove LightArtifact - default 3
entry_position+=1

SetDimLabel 0,entry_position,Brightness_cut_8bit,OS_Parameters
OS_Parameters[%Brightness_cut_8bit] = 0.01 // 0:1 - at what % brightness are 8 bit transforms cut. default 0.01 = 1%
entry_position+=1

SetDimLabel 0,entry_position,StimulatorDelay,OS_Parameters
//get some variables, 
//like the computer name (indicates which setup was used)
pcName = FindDimLabel(wParamsStr,0,"ComputerName" )

//print setup
setup = wParamsStr['pcName']
//and date of the recording
recDay = FindDimLabel(wParamsStr,0,"DateStamp_d_m_y" )
recDate = wParamsStr['recDay']
//convert each part of the string to numbers
year = str2num(recDate[0,3])
month=str2num(recDate[5,7])
day = str2num(recDate[7,9])

if (stringmatch(setup,"euler14_01")) // SETUP 1
		OS_Parameters[%StimulatorDelay] = 0// nMilliseconds delay of the stimulator between the trigger and the 
   											 // light actually hitting the tissue. for Arduino stimulator this is 0ms
   											  
elseif (stringmatch(setup,"euler14_lab2-1"))		  // SETUP 2
	if (year < 2016)
		OS_Parameters[%StimulatorDelay] = 27  // old stimulator software (using directx) written in Pascal is 27,
	elseif (year == 2016)
		if (month <= 4 && day <= 24)
			OS_Parameters[%StimulatorDelay] = 27  // old stimulator software (using directx) written in Pascal is 27,
		else
			OS_Parameters[%StimulatorDelay] = 100//new python software (using openGL) installed 25th April 2016) is 100 ms
		endif
	else
		OS_Parameters[%StimulatorDelay] = 100//new python software (using openGL) installed 25th April 2016) is 100 ms
	endif
else										  // SETUP 3 (downstairs)
	OS_Parameters[%StimulatorDelay] = 0  // Right now only has Arduino stimulator on it
endif
entry_position+=1

/// DETREND ////////////////////////////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,Detrend_smooth_window,OS_Parameters
OS_Parameters[%Detrend_smooth_window] = 1000 // smoothing window in seconds - default 1000
entry_position+=1

SetDimLabel 0,entry_position,Detrend_RatiometricData,OS_Parameters
OS_Parameters[%Detrend_RatiometricData] = 0 // Does Ratiometric data get detrended (1) or just combined (0)? - default 0
entry_position+=1

SetDimLabel 0,entry_position,Detrend_Skip,OS_Parameters
OS_Parameters[%Detrend_Skip] = 0 // 1 skips the detrend
entry_position+=1

/// MULTIPLANE ////////////////////////////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,nPlanes,OS_Parameters
OS_Parameters[%nPlanes] = 1 // number of scan planes to be deinterleaved
entry_position+=1

/// ROI PLACEMENT /////////////////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,ROI_corr_min,OS_Parameters
OS_Parameters[%ROI_corr_min] = 0.1 // Activity correlation minimum to allow a seeded ROI to grow - default 0.1
entry_position+=1

SetDimLabel 0,entry_position,ROI_mindiameter,OS_Parameters
OS_Parameters[%ROI_mindiameter] = 1 // min circle equivalent diameter in micron - default 1
entry_position+=1

SetDimLabel 0,entry_position,ROI_maxdiameter,OS_Parameters
OS_Parameters[%ROI_maxdiameter] = 3 // max circle equivalent diameter in micron - default 3
entry_position+=1

SetDimLabel 0,entry_position,ROI_minpix,OS_Parameters
OS_Parameters[%ROI_minpix] = 2 // minimum number of pixels per ROI, overrides ROI_mindiameter - default 2
entry_position+=1

SetDimLabel 0,entry_position,nRoiKillsAllowed,OS_Parameters
OS_Parameters[%nRoiKillsAllowed] = 10 // nRois that get killed due to size minimum before the routine aborts trying to place more  - default 10
entry_position+=1

SetDimLabel 0,entry_position,ROI_PxBinning,OS_Parameters
OS_Parameters[%ROI_PxBinning] = 1 // Bin pixels to autoplace ROIs (speedup = 2^Bin)  - default 1
entry_position+=1

/// TRACE AND TRIGGER EXTRACTION  ///////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,Trigger_Threshold,OS_Parameters
OS_Parameters[%Trigger_Threshold] = 20000 // Threshold to Trigger in Triggerchannel - default 20000
entry_position+=1

SetDimLabel 0,entry_position,Trigger_after_skip_s,OS_Parameters
OS_Parameters[%Trigger_after_skip_s] = 0.1 // if triggers in triggerchannel, it skips X seconds - default 0.1
entry_position+=1

SetDimLabel 0,entry_position,Trigger_LevelRead_after_lines,OS_Parameters
OS_Parameters[%Trigger_levelread_after_lines] = 2  // to read "Triggervalue" - want to avoid landing on the slope of the trigger - default 2
entry_position+=1

SetDimLabel 0,entry_position,Skip_First_Triggers,OS_Parameters 
OS_Parameters[%Skip_First_Triggers] = 0  // skips last trigger, e.g. when last loop is not complete - default 0
entry_position+=1

SetDimLabel 0,entry_position,Skip_Last_Triggers,OS_Parameters // KF 20160310
OS_Parameters[%Skip_Last_Triggers] = 0  // skips last trigger, e.g. when last loop is not complete - default 0
entry_position+=1

SetDimLabel 0,entry_position,Trigger_DisplayHeight,OS_Parameters
OS_Parameters[%Trigger_DisplayHeight] = 6  // How long are the trigger lines in the display (in SD) - default 6
entry_position+=1

SetDimLabel 0,entry_position,Baseline_nSeconds,OS_Parameters
OS_Parameters[%Baseline_nSeconds] = 5  // takes the 1st n seconds to calculate the baseline noise (for z-normalisation) - default 5
entry_position+=1

SetDimLabel 0,entry_position,Ignore1stXseconds,OS_Parameters
OS_Parameters[%Ignore1stXseconds] = 1 // for baseline extraction & for averaging across triggers (below): ignores X 1st seconds of triggers
entry_position+=1

SetDimLabel 0,entry_position,IgnoreLastXseconds,OS_Parameters
OS_Parameters[%IgnoreLastXseconds] = 0 // if weird stuff happens at end of trace can cut away
entry_position+=1

/// BASIC AVERAGING  /////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,Use_Znorm,OS_Parameters
OS_Parameters[%Use_Znorm] = 1 // use znormalised or raw traces (0/1) - default 1
entry_position+=1

SetDimLabel 0,entry_position,Trigger_Mode,OS_Parameters
OS_Parameters[%Trigger_Mode] = 1 // Use every nth trigger - default 1
entry_position+=1

SetDimLabel 0,entry_position,AverageStack_make,OS_Parameters
OS_Parameters[%AverageStack_make] = 0 // yes or no /0/1 - default 0
entry_position+=1

/// QuickCluster  /////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,Clustering_nClasses,OS_Parameters
OS_Parameters[%Clustering_nClasses] = 10 // how many Classes seeded into Clustering (kMeans)
entry_position+=1

SetDimLabel 0,entry_position,Clustering_SDplot,OS_Parameters
OS_Parameters[%Clustering_SDplot] = 5 // nSDs plotted in Clustering display
entry_position+=1


/// EVENT TRIGGERING  ////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,Events_nMax,OS_Parameters
OS_Parameters[%Events_nMax] = 1000 // maximal number of events identified in single full trace - default 1000
entry_position+=1

SetDimLabel 0,entry_position,Events_Threshold,OS_Parameters
OS_Parameters[%Events_Threshold] = 1 // Threshold for Peak detection (log scale), default = 1
entry_position+=1

SetDimLabel 0,entry_position,Events_RateBins_s,OS_Parameters
OS_Parameters[%Events_RateBins_s] = 0.05 // "Smooth_size" for Event rate plots (s) - default 0.05
entry_position+=1

/// RF Calculations  /////////////////////////////////////////////////////////////////////////

SetDimLabel 0,entry_position,Noise_PxSize_microns,OS_Parameters
OS_Parameters[%Noise_PxSize_microns] = 20 // pixel size of 3D noise - default 20 microns
entry_position+=1

SetDimLabel 0,entry_position,Noise_EventSD,OS_Parameters
OS_Parameters[%Noise_EventSD] = 0.7 // Sensitivity of Event triggering
entry_position+=1

SetDimLabel 0,entry_position,Noise_interval_sec,OS_Parameters
OS_Parameters[%Noise_interval_sec] = 0.078 // Refresh rate of the Noise (in seconds)
entry_position+=1

SetDimLabel 0,entry_position,Noise_FilterLength_s,OS_Parameters
OS_Parameters[%Noise_FilterLength_s] = 1.3 // Length extracted in seconds
entry_position+=1

SetDimLabel 0,entry_position,Kernel_SDplot,OS_Parameters
OS_Parameters[%Kernel_SDplot] = 30 // nSDs plotted in kernel function
entry_position+=1

SetDimLabel 0,entry_position,Kernel_MapPxbyPx,OS_Parameters
OS_Parameters[%Kernel_MapPxbyPx] = 1 // Kernel map pixel by pixel (1) or bigger ROIs? (0)
entry_position+=1

SetDimLabel 0,entry_position,Kernel_MapSmth,OS_Parameters
OS_Parameters[%Kernel_MapSmth] = 0 // nMicrons that kernelmaps get smoothed by
entry_position+=1

SetDimLabel 0,entry_position,Kernel_MapSDCut,OS_Parameters
OS_Parameters[%Kernel_MapSDCut] = 2 // Kills ROIs below this SD in maps
entry_position+=1

SetDimLabel 0,entry_position,Kernel_MapRange,OS_Parameters
OS_Parameters[%Kernel_MapRange] = 5 // SDs used in the RGBU maps
entry_position+=1

SetDimLabel 0,entry_position,Kernel_FFTRange,OS_Parameters
OS_Parameters[%Kernel_FFTRange] = 3 // SDs used in the FFT maps
entry_position+=1

SetDimLabel 0,entry_position,Kernel_FFTOffset,OS_Parameters
OS_Parameters[%Kernel_FFTOffset] = 1 // sets n Hz as darkest shade (default: 1Hz)
entry_position+=1

SetDimLabel 0,entry_position,Noise_Compression,OS_Parameters
OS_Parameters[%Noise_Compression] = 10 // Noise RF calculation speed up
entry_position+=1

SetDimLabel 0,entry_position,ROIKernelSmooth_space,OS_Parameters
OS_Parameters[%ROIKernelSmooth_space] = 1 // ROIKernel Smoothing (pixels)
entry_position+=1

SetDimLabel 0,entry_position,ROIKernel_seedROI,OS_Parameters
OS_Parameters[%ROIKernel_seedROI] = 0 // Which ROI is all calculated relative to, default: 0
entry_position+=1

/// redimension the OS_parameter table, so it doesn't have trailing NaN's
redimension /N=(entry_position) OS_Parameters

endif
		
// Display the Table
edit /k=1 /W=(50,50,300,700)OS_Parameters.l, OS_Parameters

end

//////////////////////////

function OS_ParameterTable_kill()

// 1 // check for Parameter Table
if (waveexists($"OS_Parameters")==1)
	wave OS_parameters
	killwaves OS_parameters	
	print "OS_Parameters deleted"
else
	print "OS_Parameters was already deleted"
endif
end
