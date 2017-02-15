#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_LineScanFormat()

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
variable X_cut = OS_Parameters[%LightArtifact_cut]
variable LineDuration = OS_Parameters[%LineDuration]
variable Display_LineImage = OS_Parameters[%Display_Stuff]

printf "Converting wDataCh"
printf Num2Str(Channel)
printf " to linescan stack..."

// data handling
string input_name = "wDataCh"+Num2Str(Channel)
string input_name2 = "wDataCh"+Num2Str(TriggerChannel)
duplicate /o $input_name InputData
duplicate /o $input_name2 InputTriggers


variable nX = DimSize(InputData,0)
variable nY = DimSize(InputData,1)
variable nF = DimSize(InputData,2)

if (nY>1)
	duplicate /o $input_name RawData_original
	duplicate /o $input_name2 RawTriggers_original
	print "The original wDataCh waves are saved as RawData_original & RawTriggers_original"
endif

// calculate Pixel size in microns to scale ROIs
wave wParamsNum // Reads data-header
variable zoom = wParamsNum(30) // extract zoom
variable px_Size = (0.65/zoom * 110)/nX // microns

// make stack with single line in Y
make /o/n=(nX,1,nY*nF) OutputData = NaN
make /o/n=(nX,1,nY*nF) OutputTriggers = NaN
make /o/n=(nY*nF,nX-X_cut) LineScanImage = NaN
setscale /p y,-nX/2*px_Size,px_Size,"µm" LineScanImage
Setscale /p x,0,LineDuration,"s" LineScanImage


variable yy,ff
for (ff=0; ff<nF; ff+=1)
	for (yy=0; yy<nY; yy+=1)
		Multithread OutputData[][0][ff*nY+yy]=InputData[p][yy][ff]
		Multithread OutputTriggers[][0][ff*nY+yy]=InputTriggers[p][yy][ff]
		Multithread LineScanImage[ff*nY+yy][]=InputData[q+X_cut][yy][ff]
	endfor
endfor
duplicate /o OutputData $input_name // overwrites the wDataChX wave
duplicate /o OutputTriggers $input_name2 // overwrites the wDataChX wave
killwaves OutputData,InputData,InputTriggers
print " complete..."


// Display

if (Display_LineImage==1)
	Display /k=1
	Appendimage LineScanIMage
	ModifyGraph fSize=8,axisEnab(left)={0.05,1},axisEnab(bottom)={0.05,1};DelayUpdate
	Label bottom "\\Z10Time (\\U)"
endif

end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////