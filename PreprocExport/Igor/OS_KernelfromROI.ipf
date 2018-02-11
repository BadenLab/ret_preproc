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
variable SeedROI = OS_Parameters[%ROIKernel_seedROI]

// hardcoded for now: - should eventually go into OS_Parameters
variable nS_Montage_pre = 0.5 // s - specifically for Montage - the full kernel dur is set by FilterHalfLength_s above
variable nS_Montage_post = 1 // s
variable nS_CrossCorr = 0.4 // s - total duration, i.e. half prior half after event
variable CrossCorr_UseDiF = 0 // does crosscorr use Differentiated traces (1) or raw (0)?
variable DiffSmooth = 20 // now many points (lines) are crosscrorr input traces smoothed
variable ROIKernelSmooth_time = 2

// data handling
wave wParamsNum // Reads data-header
wave ROIs
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
string output_name4 = "CrossCorrStack"+Num2Str(Channel)
string output_name5 = "ROIKernel_traces"+Num2Str(Channel)
string output_name6 = "ROICrossCorr_traces"+Num2Str(Channel)

make /o/n=(nX, nY) ROIs_XCrop = ROIs[p+X_Cut][q]


variable xx,yy,ff,rr

///////////////////// MAIN //////////////////////////////////////////////////////

variable Skip1stNFrames = Ignore1stXseconds * Framerate
variable SkipLastNFrames = IgnoreLastXseconds * Framerate

variable nF_KernelNeg = FilterHalfLength_s * Framerate
variable nF_KernelPos = FilterHalfLength_s * Framerate
variable nFrames_baseline = nF_KernelNeg/2
variable nF_Kernel = nF_KernelNeg + nF_KernelPos
variable skipframes = 1

make /o/n=(nF) SeedTrace = InputTraces[p][SeedROI]
Differentiate/DIM=0  SeedTrace/D=SeedTrace_DIF
WaveStats/Q SeedTrace_Dif
SeedTrace_Dif/=V_SDev

variable nF_Montage_pre = ceil(nS_Montage_pre * Framerate)
variable nF_Montage_post = ceil(nS_Montage_post * Framerate)
variable nF_Montage = nF_Montage_pre+nF_Montage_post




variable nP_CrossCorr = (nS_CrossCorr / LineDuration) 

make /o/n=(nX,nY,nF_Kernel) ROIKernelStack = 0
make /o/n=(nX,nY,nP_CrossCorr) CrossCorrStack = 0
make /o/n=(nX,nY) ROIKernelImage = 0
make /o/n=(nX,(nY+2)*nF_Montage+2) ROIKernelMontage = 0

setscale /p z,-nS_CrossCorr/2,LineDuration,"s" CrossCorrStack


// generate kernelstack
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

// generate CrossCorr stack

Smooth DiffSmooth, SeedTrace_DIF
for (xx=0;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		make /o/n=(nF) TargetTrace = InputStack[xx][yy][p]
		Differentiate/DIM=0  TargetTrace/D=TargetTrace_DIF
		Smooth DiffSmooth, TargetTrace_DIF
		Correlate/NODC SeedTrace_DIF, TargetTrace_DIF
		Correlate/NODC SeedTrace, TargetTrace
		if (CrossCorr_UseDiF==1)
			Multithread CrossCorrStack[xx][yy][]=TargetTrace_DIF[r+nF -nP_CrossCorr/2]
		else
			Multithread CrossCorrStack[xx][yy][]=TargetTrace_DIF[r+nF -nP_CrossCorr/2]
		endif
	endfor
endfor

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
if (ROIKernelSmooth_time>0)
	Smooth /DIM=2 ROIKernelSmooth_time, ROIKernelStack
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
	ROIKernelMontage[][ff*(nY+2),(ff+1)*(nY+2)-3][ff]=ROIKernelStack[p][q-ff*(nY+2)][ff+nF_KernelNeg-nF_Montage_pre]
endfor

// Get Kernel and CrossCorr Traces from ROIs
make /o/n=(nF_Kernel,nROIs) ROIKernel_Traces = 0
make /o/n=(nP_CrossCorr,nROIs) ROICrossCorr_Traces = 0

for (rr=0;rr<nRois;rr+=1)
	variable ROI_value = (rr+1)*-1 // ROIs in Mask are coded as negative starting from -1 (SARFIA standard)
	variable ROI_size = 0
	for (xx=0;xx<nX;xx+=1)
		for (yy=0;yy<nY;yy+=1)
			if (ROIs[xx][yy]==ROI_value)
				ROI_size+=1
				ROIKernel_Traces[][rr]+=ROIKernelStack[xx][yy][p] // add up each pixel of a ROI
				ROICrossCorr_Traces[][rr]+=CrossCorrStack[xx][yy][p] // add up each pixel of a ROI
			endif
		endfor
	endfor
	ROIKernel_Traces[][rr]/=ROI_size // now is average activity of ROI
	ROICrossCorr_Traces[][rr]/=ROI_size // now is average activity of ROI	

	make /o/n=(nP_CrossCorr) currentwave = ROICrossCorr_Traces[p][rr]
	WaveStats/Q currentwave
	ROICrossCorr_Traces[][rr]/=V_Max
	
	make /o/n=(nF_Kernel) currentwave = ROIKernel_Traces[p][rr]
	WaveStats/Q currentwave
	ROIKernel_Traces[][rr]/=V_Max

endfor



//////////////
setscale /p x,-nX/2*px_Size,px_Size,"µm" ROIKernelImage, ROIKernelStack, ROIKernelMontage, ROIs_XCrop
setscale /p y,-nY/2*px_Size,px_Size,"µm"  ROIKernelImage, ROIKernelStack, ROIs_XCrop
setscale y,-nF_Montage_pre*nY*LineDuration,nF_Montage_post*nY*LineDuration,"s"  ROIKernelMontage
setscale /p x,-nS_CrossCorr/2,LineDuration,"s" ROICrossCorr_Traces
setscale /p x,-FilterHalfLength_s,LineDuration * nY,"s" ROIKernel_Traces


// export handling
duplicate /o ROIKernelStack $output_name1
duplicate /o ROIKernelImage $output_name2
duplicate /o ROIKernelMontage $output_name3
duplicate /o CrossCorrStack $output_name4
duplicate /o ROIKernel_traces $output_name5
duplicate /o ROICrossCorr_traces $output_name6

	
// display

if (Display_stuff==1)

	// MONTAGE

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

	Appendimage /l=image2Y /b=image2X ROIs_XCrop
	ModifyImage ROIs_XCrop explicit=1,eval={-SeedROI-1,65280,65280,65280}

	ModifyGraph swapXY=1
	DoUpdate
	ModifyGraph width=0,height=0
	
	// Traces
	Display /k=1
	ModifyGraph height={Aspect,2.5}
	ModifyGraph width=150
	Appendimage /l=imageY /b=imageX $output_name2
	Appendimage /l=imageY /b=imageX ROIs_XCrop
	
		// colour in the ROIs
		make /o/n=(1) M_Colors
		Colortab2Wave Rainbow256
		for (rr=0;rr<nRois;rr+=1)
			variable colorposition = 255 * (rr+1)/nRois
			ModifyImage ROIs_XCrop explicit=1,eval={-rr-1,M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2]}
		endfor
		ModifyImage ROIs_XCrop explicit=1,eval={-SeedROI-1,65280,65280,65280}
	
	
	
	ModifyGraph noLabel(imageY)=2,noLabel(imageX)=2,axThick(imageY)=0;DelayUpdate
	ModifyGraph axThick(imageX)=0,axisEnab(imageY)={0.85,1},axisEnab(imageX)={0.05,1};DelayUpdate
	ModifyGraph freePos(imageY)={0,kwFraction},freePos(imageX)={0.8,kwFraction}
	
	for (rr=0;rr<nROIs;rr+=1)
		string tracename1 = output_name5+"#"+Num2Str(rr)
		string tracename2 = output_name6+"#"+Num2Str(rr)		
		if (rr==0)
			tracename1 = output_name5
			tracename2 = output_name6
		endif
		Appendtograph /l=KernelTraceY /b=KernelTraceX $output_name5[][rr]
		Appendtograph /l=CrossCorrTraceY /b=CrossCorrX $output_name6[][rr]

		colorposition = 255 * (rr+1)/nRois
		ModifyGraph rgb($tracename1)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])	
		ModifyGraph rgb($tracename2)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])			
	
	endfor
	
		// kill the seed traces
		tracename1 = output_name5+"#"+Num2Str(SeedROI)
		tracename2 = output_name6+"#"+Num2Str(SeedROI)		
		if (SeedROI==0)
			tracename1 = output_name5
			tracename2 = output_name6
		endif
		RemoveFromGraph $tracename1
		RemoveFromGraph $tracename2
	
	
	ModifyGraph fSize=8,lblPos=47,axisEnab(KernelTraceY)={0.05,0.3};DelayUpdate
	ModifyGraph axisEnab(KernelTraceX)={0.05,1},axisEnab(CrossCorrTraceY)={0.5,0.75};DelayUpdate
	ModifyGraph axisEnab(CrossCorrX)={0.05,1},freePos(KernelTraceY)={0,kwFraction};DelayUpdate
	ModifyGraph freePos(KernelTraceX)={0,kwFraction};DelayUpdate
	ModifyGraph freePos(CrossCorrTraceY)={0,kwFraction};DelayUpdate
	ModifyGraph freePos(CrossCorrX)={0.45,kwFraction};DelayUpdate
	Label KernelTraceY "\\Z10ROI Kernel (z-scores)";DelayUpdate
	Label KernelTraceX "\\Z10Time (\\U)";DelayUpdate
	Label CrossCorrTraceY "\\Z10CrossCorr (A.U.)";DelayUpdate
	Label CrossCorrX "\\Z10Time (\\U)"
	ModifyGraph lblPos(CrossCorrX)=37
	ModifyGraph lsize=1.5
	
	DoUpdate
	ModifyGraph width=0,height=0
	
	
	
	
endif


// cleanup
killwaves currentwave, ROIKernelStack, ROIKernelImage, InputStack, SeedTrace_DIF, SeedTrace, InputTraces, ROIKernelMontage
killwaves CrossCorrStack, ROIKernel_traces, ROICrossCorr_traces, TargetTrace_DIF, TargetTrace




end