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
variable nSeconds_smooth = OS_Parameters[%Detrend_smooth_window]
variable LightArtifactCut = OS_Parameters[%LightArtifact_cut]

// data handling
string input_name = "wDataCh"+Num2Str(Channel)
string output_name = "wDataCh"+Num2Str(Channel)+"_detrended"
duplicate /o $input_name InputData
variable nX = DimSize(InputData,0)
variable nY = DimSize(InputData,1)
variable nF = DimSize(InputData,2)
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
Smooth/DIM=2 Smoothingfactor, InputData
Multithread OutputData[][][]-=InputData[p][q][r]-Mean_image[p][q]

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
killwaves CurrentTrace,InputData,OutputData,mean_image

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

