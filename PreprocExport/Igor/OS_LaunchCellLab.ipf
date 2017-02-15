#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "CellLabUI"
#include "CellLabRoutines"

/////////////////////////////////////////////////////////////////////////////////////////////////////
///	Official ScanM Data Preprocessing Scripts - by Tom Baden    	///
/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_LaunchCellLab()

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
variable Display_RoiMask = OS_Parameters[%Display_Stuff]
variable X_cut = OS_Parameters[%LightArtifact_cut]

// data handling
string input_name = "wDataCh"+Num2Str(Channel)+"_detrended"
duplicate /o $input_name InputData
variable nX = DimSize(InputData,0)
variable nY = DimSize(InputData,1)
variable nF = DimSize(InputData,2)
//variable Framerate = 1/(nY * LineDuration) // Hz 
//variable Total_time = (nF * nX ) * LineDuration
//print "Recorded ", total_time, "s @", framerate, "Hz"
variable xx,yy

// make average image
make /o/n=(nX,nY) Stack_ave = InputData[nX/2][nY/2][nF/2] // so that Light Artifact is rioughly the same brightness as the rest of the scan
for (xx=X_Cut;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		make /o/n=(nF) currentwave = InputData[xx][yy][p]
		Wavestats/Q Currentwave
		Stack_ave[xx][yy]=V_Avg
	endfor
endfor




// CellLab only accepts square scans. So if scan is not square, here extend the Image it works with to be square
if (nX==nY)
	string sourcename = "Stack_ave"

else
	make /o/n=(nX,nX) Stack_Ave_square = InputData[nX/2][nY/2][nF/2]
	Stack_Ave_square[][0,nX-nY]=Stack_Ave[p][q]
	//Setscale 
	sourcename = "Stack_ave_square"
	
endif
	string targetname = "ROIs"

// call cell lab

RecognizeCellsUI($sourcename, targetname, "")




// cleanup
killwaves currentwave, InputData

end