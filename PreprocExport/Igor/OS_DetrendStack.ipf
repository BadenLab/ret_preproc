#pragma rtGlobals=3		// Use modern global access method and strict wave access.

/////////////////////////////////////////////////////////////////////////////////////////////////////
///	Official ScanM Data Preprocessing Scripts - by Tom Baden    	///
/////////////////////////////////////////////////////////////////////////////////////////////////////
///	Requires raw 3D data in 16 bits with no preprocessing			///
///	Input Arguments - which Channel (0,1,2...?)				     	///
///	e.g. "OS_DetrendStack(0)	"							     	///
///   --> reads wDataCh0,1,2...									///
///   --> for each pixel subtracts heavily smoothed version of itself   	///
///   --> ...and adds its own mean (to avoid going out of range)		///
///	Output is new wave called wDataCh..._detrended				///
/////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_DetrendStack()

printf "Detrending...."

// flags from "OS_Parameters"
if (waveexists($"OS_Parameters")==0)
	print "Warning: OS_Parameters wave not yet generated - doing that now..."
	OS_ParameterTable()
	DoUpdate
endif
wave OS_Parameters
variable Channel = OS_Parameters[%Data_Channel]
variable TriggerChannel = OS_Parameters[%Trigger_Channel]
variable nSeconds_smooth = OS_Parameters[%Detrend_smooth_window]
variable LightArtifactCut = OS_Parameters[%LightArtifact_cut]
variable nPlanes = OS_Parameters[%nPlanes]
variable SkipDetrend = OS_Parameters[%Detrend_skip]

// data handling
string input_name = "wDataCh"+Num2Str(Channel)
string input_name2 = "wDataCh"+Num2Str(TriggerChannel)
string output_name = "wDataCh"+Num2Str(Channel)+"_detrended"
string output_name2 = "wDataCh"+Num2Str(TriggerChannel) // this will overwrite wDataCh2 e.g. - but there is the TriggerData_original backup

duplicate /o $input_name InputData_preflip
if (waveexists($"TriggerData_original")==0)
	duplicate /o $input_name2 TriggerData_original
	duplicate /o $input_name2 TriggerData
else
	wave TriggerData_original
	duplicate /o TriggerData_original TriggerData
endif



variable nX = DimSize(InputData_preflip,0)
variable nY = DimSize(InputData_preflip,1)
variable nF = DimSize(InputData_preflip,2)
variable pp

// X-Flip InputStack, but spare the light Artifact
duplicate /o InputData_preflip InputData
InputData[LightArtifactCut,nX-1][][]=InputData_preflip[nX-1-(p-LightArtifactCut)][q][r]
killwaves InputData_preflip

// multiplane deinterleave
if (nPlanes>1)
	print "Deinterleaving", nPlanes, "planes"
	
	variable nF_true = floor(nF / nPlanes)
	
	make /o/n=(nX, nY*nPlanes, nF_true) InputData_deinterleaved = NaN
	make /o/n=(nX, nY*nPlanes, nF_true) TriggerData_deinterleaved = NaN
	for (pp=0;pp<nPlanes;pp+=1)
		Multithread InputData_deinterleaved[][nY*pp,nY*(pp+1)-1][]=InputData[p][q-nY*pp][r*nPlanes+pp]
		Multithread TriggerData_deinterleaved[][nY*pp,nY*(pp+1)-1][]=TriggerData[p][q-nY*pp][r*nPlanes+pp]
	endfor
	Duplicate/o InputData_deinterleaved InputData
	Duplicate/o TriggerData_deinterleaved TriggerData	
	nY*=nPlanes
	nF=nF_true

endif
duplicate/o InputData OutputData

// calculate size of smoothing window
variable Framerate = 1/(nY * 0.002) // Hz
variable Smoothingfactor = Framerate * nSeconds_smooth
if (Smoothingfactor>2^15-1) // exception handling - limit smooth function to its largest allowed input
	Smoothingfactor = 2^15-1 
endif

// detrending
variable xx,yy
make /o/n=(nX,nY) mean_image = 0
for (xx=0; xx<nX; xx+=1)
	for (yy=0; yy<nY; yy+=1)
		make/o/n=(nF) CurrentTrace = InputData[xx][yy][p]
		Wavestats/Q CurrentTrace
		mean_image[xx][yy]=V_Avg
	endfor
endfor
if (SkipDetrend==0)
	Smooth/DIM=2 Smoothingfactor, InputData
	Multithread OutputData[][][]-=InputData[p][q][r]-Mean_image[p][q]
else
	Multithread OutputData[][][]=InputData[p][q][r]
endif



// cut things
OutputData[][][0]=OutputData[p][q][1] // copy second frame into 1st to kill frame 1 artifact
make /o/n=(nX-LightArtifactCut,nY) tempimage = Mean_image[p+LightArtifactCut][q]
ImageStats/Q tempimage
OutputData[0,LightArtifactCut][][]=V_Avg // Clip Light Artifact
Mean_image[0,LightArtifactCut][]=V_Avg
killwaves tempimage

// generate output
duplicate /o OutputData $output_name
duplicate /o TriggerData $output_name2
duplicate /o mean_image Stack_Ave

// cleanup
killwaves CurrentTrace,InputData,OutputData,mean_image, InputData_deinterleaved, TriggerData_deinterleaved, TriggerData

// outgoing dialogue
print " complete..."

end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_DetrendRatiometric()

// flags from "OS_Parameters"
if (waveexists($"OS_Parameters")==0)
	print "Warning: OS_Parameters wave not yet generated - doing that now..."
	OS_ParameterTable()
	DoUpdate
endif
wave OS_Parameters
variable Channel = OS_Parameters[%Data_Channel]
variable Channel2 = OS_Parameters[%Data_Channel2]
variable nSeconds_smooth = OS_Parameters[%Detrend_smooth_window]
variable Detrend_Ratiometricdata = OS_Parameters[%Detrend_RatioMetricData]
variable LightArtifactCut = OS_Parameters[%LightArtifact_cut]

// data handling
string input_name = "wDataCh"+Num2Str(Channel)
string input_name2 = "wDataCh"+Num2Str(Channel2)
string output_name = "wDataCh"+Num2Str(Channel)+"_detrended"
duplicate /o $input_name InputData
duplicate /o $input_name2 InputData2
variable nX = DimSize(InputData,0)
variable nY = DimSize(InputData,1)
variable nF = DimSize(InputData,2)
duplicate/o InputData OutputData

// Get RatioMetric Stack (after which everything is identical to DetrendStack routine)
make /o/n=(nX,nY) InputData2_frame2 = InputData2[p][q][1]
ImageStats/Q InputData2_frame2
variable InputData2_brightness = V_Avg
InputData[][][]/=InputData2[p][q][r]/InputData2_brightness

if (Detrend_RatiometricData==0)
	print "Complete... (no detrending done, only Channel division)"
	OutputData[][][]=InputData[p][q][r]
else
	// calculate size of smoothing window
	variable Framerate = 1/(nY * 0.002) // Hz
	variable Smoothingfactor = Framerate * nSeconds_smooth
	if (Smoothingfactor>2^15-1) // exception handling - limit smooth function to its largest allowed input
		Smoothingfactor = 2^15-1 
	endif
	
	// detrending
	printf "Detrending...."
	variable xx,yy
	make /o/n=(nX,nY) mean_image = 0
	for (xx=0; xx<nX; xx+=1)
		for (yy=0; yy<nY; yy+=1)
			make/o/n=(nF) CurrentTrace = InputData[xx][yy][p]
			Wavestats/Q CurrentTrace
			mean_image[xx][yy]=V_Avg
		endfor
	endfor
	Smooth/DIM=2 Smoothingfactor, InputData
	Multithread OutputData[][][]-=InputData[p][q][r]-Mean_image[p][q]
	print " complete..."
endif	

// cut things
OutputData[][][0]=OutputData[p][q][1] // copy second frame into 1st to kill frame 1 artifact
make /o/n=(nX-LightArtifactCut,nY) tempimage = Mean_image[p+LightArtifactCut][q]
ImageStats/Q tempimage
OutputData[0,LightArtifactCut][][]=V_Avg // Clip Light Artifact
Mean_image[0,LightArtifactCut][]=V_Avg
killwaves tempimage
	
// generate output
duplicate /o OutputData $output_name
duplicate /o mean_image Stack_Ave

// cleanup
killwaves CurrentTrace,InputData,OutputData,InputData2,InputData2_frame2,Mean_image
	

	
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_SaveRawAsTiff()

printf "Saving PreProc Stack as Tiff..."

// 1 // check for Parameter Table
if (waveexists($"OS_Parameters")==0)
	print "Warning: OS_Parameters wave not yet generated - doing that now..."
	OS_ParameterTable()
	DoUpdate
endif
wave OS_Parameters
// 2 //  check for Detrended Data stack
variable Channel = OS_Parameters[%Data_Channel]
if (waveexists($"wDataCh"+Num2Str(Channel)+"_detrended")==0)
	print "Warning: wDataCh"+Num2Str(Channel)+"_detrended wave not yet generated - doing that now..."
	OS_DetrendStack()
endif

// flags from "OS_Parameters"
variable X_cut = OS_Parameters[%LightArtifact_cut]
variable brightness_cutoff = OS_Parameters[%Brightness_cut_8bit]

// data handling
string input_name = "wDataCh"+Num2Str(Channel)+"_detrended"
duplicate /o $input_name InputData
variable nX = DimSize(InputData,0) - X_cut
variable nY = DimSize(InputData,1)
variable nF = DimSize(InputData,2)

string output_name = "wDataCh"+Num2Str(Channel)+"_detrended_8bit"

variable ff, xx, yy,pp

// Convert to 8 bit & flip Y axis (as imageJ is default Y flipped relative to Igor)
make /o/n=(nX,nY,nF) OutputData = InputData[p+X_Cut][nY-1-q][r]

	// get brightness histogram
Make/N=(2^16) /O Brightness_Hist = 0
Make/N=(2^16) /O currentwave_Hist = 0
make /o/n=(nF) currentwave = 0
for (xx=0;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		Multithread currentwave = InputData[xx+X_cut][yy][p]
		Histogram/B={0,1,2^16} currentwave,currentwave_Hist
		Multithread Brightness_Hist[]+=currentwave_Hist[p]
	endfor
endfor
WaveStats/Q Brightness_Hist
Brightness_Hist/=V_Max

	// find start of brightness distribution
variable Min_brightness = 0
for (pp=0; pp<2^16;pp+=1)
	if (Brightness_Hist[pp]>brightness_cutoff)
		Min_brightness = pp
		pp = 2^16 // abort
	endif
endfor
	// find end of brightness distribution
variable Max_brightness = Min_brightness
for (pp=Min_brightness; pp<2^16;pp+=1)
	if (Brightness_Hist[pp]<brightness_cutoff)
		Max_brightness = pp
		pp = 2^16 // abort
	endif
endfor
print "Setting brightness between", Min_brightness, "and", Max_brightness

// scaling outputwave
OutputData[][][]=(OutputData[p][q][r]<Min_brightness)?(Min_brightness):(OutputData[p][q][r])
OutputData[][][]=(OutputData[p][q][r]>Max_brightness)?(Max_brightness):(OutputData[p][q][r])

OutputData[][][]-=Min_brightness
OutputData[][][]/=(Max_brightness-Min_brightness) / (2^8)




// generate output
duplicate /o OutputData $output_name
imagesave /s/t="tiff" $output_name


// cleanup
killwaves inputData, OutputData,currentwave,currentwave_Hist

// outgoing dialogue
print " complete..."

end