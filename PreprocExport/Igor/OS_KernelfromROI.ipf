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
variable Kernel_MapRange = OS_Parameters[%Kernel_MapRange]

variable nS_Montage_pre = 0.5
variable nS_Montage_post = 1


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
string output_name3 = "ROIKernelMontage"+Num2Str(Channel)


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

variable nF_Montage_pre = ceil(nS_Montage_pre * Framerate)
variable nF_Montage_post = ceil(nS_Montage_post * Framerate)
variable nF_Montage = nF_Montage_pre+nF_Montage_post

make /o/n=(nX,nY,nF_Kernel) ROIKernelStack = 0
make /o/n=(nX,nY) ROIKernelImage = 0
make /o/n=(nX,(nY+2)*nF_Montage) ROIKernelMontage = 0

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
ROIKernelStack[][][]=(NumType(ROIKernelStack[p][q][r])==2)?(0):(ROIKernelStack[p][q][r])


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
		ROIKernelImage[xx][yy]=V_SDev
	endfor
endfor

// make Montage
for (ff=0;ff<nF_Montage;ff+=1)
	ROIKernelMontage[][ff*(nY+2),(ff+1)*(nY+2)-3][ff]=ROIKernelStack[p][q-ff*(nY+2)][ff+nF_KernelNeg-nF_Montage_pre]
endfor

//////////////
setscale /p x,-nX/2*px_Size,px_Size,"µm" ROIKernelImage, ROIKernelStack, ROIKernelMontage
setscale /p y,-nY/2*px_Size,px_Size,"µm"  ROIKernelImage, ROIKernelStack

setscale y,-nF_Montage_pre*nY*LineDuration,nF_Montage_post*nY*LineDuration,"s"  ROIKernelMontage



// export handling
duplicate /o ROIKernelStack $output_name1
duplicate /o ROIKernelImage $output_name2
duplicate /o ROIKernelMontage $output_name3

	
// display

if (Display_stuff==1)

	display /k=1
	variable Aspectratio = (nY / nX) * nF_Montage
	ModifyGraph height={Aspect,Aspectratio}
	ModifyGraph width=80
	Appendimage /l=imageY /b=imageX $output_name3
	ModifyGraph fSize=8,lblPos=47,axisEnab(imageY)={0.05,1},axisEnab(imageX)={0.05,1};DelayUpdate
	ModifyGraph freePos(imageY)={0,kwFraction},freePos(imageX)={0,kwFraction};DelayUpdate
	Label imageY "\\Z10Time (\\U)"
	Label imageX "\\Z10\\U"
	ModifyGraph zero(imageY)=1
	ModifyImage $output_name3 ctab= {0,Kernel_MapRange,VioletOrangeYellow,0}
	
	Appendimage /l=image2Y /b=image2X $output_name2
	ModifyGraph fSize(image2Y)=8,noLabel(image2Y)=2,noLabel(image2X)=2;DelayUpdate
	ModifyGraph axThick(image2Y)=0,axThick(image2X)=0,axisEnab(image2Y)={0.02,0.035};DelayUpdate
	ModifyGraph axisEnab(image2X)={0.05,1},freePos(image2Y)={0,kwFraction};DelayUpdate
	ModifyGraph freePos(image2X)={0,kwFraction}
	ModifyImage $output_name2 ctab= {0,*,VioletOrangeYellow,0}

	ModifyGraph swapXY=1
endif


// cleanup
killwaves currentwave, ROIKernelStack, ROIKernelImage, InputStack, SeedTrace_DIF, SeedTrace, InputTraces, ROIKernelMontage




end