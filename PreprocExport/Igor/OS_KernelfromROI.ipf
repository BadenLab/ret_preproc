#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_KernelFromROI()

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
// 3 //  check for ROI_Mask
if (waveexists($"ROIs")==0)
	print "Warning: ROIs wave not yet generated - doing that now (using correlation algorithm)..."
	OS_AutoRoiByCorr()
	DoUpdate
endif
// 4 //  check if Traces and Triggers are there
if (waveexists($"Triggertimes")==0)
	print "Warning: Traces and Trigger waves not yet generated - doing that now..."
	OS_TracesAndTriggers()
	DoUpdate
endif

// flags from "OS_Parameters"
variable Display_stuff = OS_Parameters[%Display_Stuff]
variable use_znorm = OS_Parameters[%Use_Znorm]
variable LineDuration = OS_Parameters[%LineDuration]
variable Ignore1stXseconds = OS_Parameters[%Ignore1stXseconds]
variable IgnoreLastXseconds = OS_Parameters[%IgnoreLastXseconds]
variable SD_threshold = OS_Parameters[%Noise_EventSD]
variable FilterHalfLength_s = OS_Parameters[%Noise_FilterLength_s]
variable X_cut = OS_Parameters[%LightArtifact_cut]
variable ROIKernelSmooth_space = OS_Parameters[%ROIKernelSmooth_space]

// data handling
wave wParamsNum // Reads data-header
string input_name = "wDataCh"+Num2Str(Channel)+"_detrended"
string traces_name = "Traces"+Num2Str(Channel)+"_raw"
if (use_znorm==1)
	traces_name = "Traces"+Num2Str(Channel)+"_znorm"
endif
duplicate /o $input_name InputStack
duplicate /o $traces_name InputTraces

variable nF = DimSize(InputTraces,0)
variable nRois = DimSize(InputTraces,1)
variable nX = DimSize(InputStack,0) - X_cut
variable nY = DimSize(InputStack,1)
variable Framerate = 1 / (nY * LineDuration)

variable zoom = wParamsNum(30) // extract zoom
variable px_Size = (0.65/zoom * 110)/nX // microns

string output_name1 = "ROIKernelStack"+Num2Str(Channel)
string output_name2 = "ROIKernelImage"+Num2Str(Channel)


variable xx,yy,ff

///////////////////// MAIN //////////////////////////////////////////////////////

variable Skip1stNFrames = Ignore1stXseconds * Framerate
variable SkipLastNFrames = IgnoreLastXseconds * Framerate

variable nF_KernelNeg = FilterHalfLength_s * Framerate
variable nF_KernelPos = FilterHalfLength_s * Framerate
variable nFrames_baseline = nF_KernelNeg/2
variable nF_Kernel = nF_KernelNeg + nF_KernelPos
variable skipframes = 1

make /o/n=(nF) SeedTrace = InputTraces[p][0]
Differentiate/DIM=0  SeedTrace/D=SeedTrace_DIF
WaveStats/Q SeedTrace_Dif
SeedTrace_Dif/=V_SDev

make /o/n=(nX,nY,nF_Kernel) ROIKernelStack = 0
make /o/n=(nX,nY) ROIKernelImage = 0

variable nTrigs = 0

for (ff=nF_KernelNeg+Skip1stNFrames;ff<nF-nF_KernelPos-SkipLastNFrames;ff+=1)
	if (SeedTrace_DIF[ff]>SD_threshold)
		Multithread ROIKernelStack[][][]+=InputStack[p+X_cut][q][ff-nF_KernelNeg+r] * (SeedTrace_DIF[ff] - SD_threshold)
		nTrigs+=1
		ff+=skipframes
	endif
endfor

ROIKernelStack/=nTrigs
Print nTrigs, "events triggered"

// znorm to baseline
make /o/n=(nFrames_baseline) currentwave = 0
for (xx=0;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		currentwave[]=ROIKernelStack[xx][yy][p]
		WaveStats/Q currentwave
		ROIKernelStack[xx][yy][]-=V_Avg
		ROIKernelStack[xx][yy][]/=V_SDev
	endfor
endfor

// make SD projection image
make /o/n=(nF_Kernel) currentwave = 0

if (ROIKernelSmooth_space>0)
	Smooth /DIM=0 ROIKernelSmooth_space, ROIKernelStack
	Smooth /DIM=1 ROIKernelSmooth_space, ROIKernelStack
endif

for (xx=0;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		currentwave[]=ROIKernelStack[xx][yy][p]
		WaveStats/Q currentwave
		ROIKernelImage[xx][yy][]=V_SDev
	endfor
endfor
//////////////
setscale /p x,-nX/2*px_Size,px_Size,"µm" ROIKernelImage, ROIKernelStack
setscale /p y,-nY/2*px_Size,px_Size,"µm"  ROIKernelImage, ROIKernelStack

// export handling
duplicate /o ROIKernelStack $output_name1
duplicate /o ROIKernelImage $output_name2

	
// display

if (Display_stuff==1)
	display /k=1
	Appendimage /l=imageY /b=imageX $output_name2
	ModifyGraph fSize=8,lblPos=47,axisEnab(imageY)={0.05,1},axisEnab(imageX)={0.05,1};DelayUpdate
	ModifyGraph freePos(imageY)={0,kwFraction},freePos(imageX)={0,kwFraction};DelayUpdate
	Label imageY "\\Z10\\U";DelayUpdate
	Label imageX "\\Z10\\U"
endif


// cleanup
killwaves currentwave, ROIKernelStack, ROIKernelImage, InputStack, SeedTrace_DIF, SeedTrace, InputTraces



end