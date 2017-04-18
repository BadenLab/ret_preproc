#pragma rtGlobals=3     // Use modern global access method and strict wave access.
function OS_LED_Noise()
// 1 // check for Parameter Table
if (waveexists($"NoiseArray4LEDs")==0)
    print "Warning: NoiseArray4LEDs wave missing - please import! Procedure aborted."
    abort
endif
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
variable Display_kernels = OS_Parameters[%Display_Stuff]
variable use_znorm = OS_Parameters[%Use_Znorm]
variable LineDuration = OS_Parameters[%LineDuration]
variable noise_interval = OS_Parameters[%Noise_interval_sec] // refresh time of Noise instances
variable Noise_Threshold = OS_Parameters[%Noise_EventSD] // nSD over baseline in time differential //read from OS_Parameters
variable nSeconds_kernel = OS_Parameters[%Noise_FilterLength_s] // nSD over baseline in time differential //read from OS_Parameters
variable nSDplot = OS_Parameters[%Kernel_SDplot] // nSD plotted in overview on y axis

// data handling
wave NoiseArray4LEDs // official stimulus array
string traces_name = "Traces"+Num2Str(Channel)+"_raw"
if (use_znorm==1)
    traces_name = "Traces"+Num2Str(Channel)+"_znorm"
endif
string tracetimes_name = "Tracetimes"+Num2Str(Channel)
duplicate /o $traces_name InputTraces
duplicate /o $tracetimes_name InputTraceTimes
wave Triggertimes
variable nF = DimSize(InputTraces,0)
variable nRois = DimSize(InputTraces,1)
string output_name1 = "Kernels"+Num2Str(Channel)
variable pp,ll,tt,rr,kk
variable nSeconds_kernel_prezero = nSeconds_kernel-0.3
variable nSeconds_kernel_baseline = 0.2
variable nSeconds_kernel_eventline = 0.8 // last X s
variable highlightSD = 2
variable suppressSD = 1
/// 
// calculating basic parameters
variable nP_stim = Dimsize(NoiseArray4LEDs,0)
variable nP_data = Dimsize(InputTraceTimes,0)
variable nLEDs = Dimsize(NoiseArray4LEDs,1)
variable nTriggers = Dimsize(Triggertimes,0)
variable timebase_s_stim =noise_interval
variable timebase_s_data = InputTraceTimes[1][0]-InputTraceTimes[0][0]
variable nP_data_upsampled = ceil(nP_data * timebase_s_data * 1/LineDuration)
variable nP_stim_upsampled = ceil(nP_stim * timebase_s_stim * 1/LineDuration)
variable nStim_repeats = ceil(nP_data_upsampled / nP_stim_upsampled )
make /o/n=(nP_stim*nStim_repeats,nLEDs) Stimulus = NaN
for (rr=0;rr<nStim_repeats;rr+=1)
    Stimulus[nP_stim*rr,nP_stim*(rr+1)-1][0,nLEDs-2]=NoiseArray4LEDs[p-nP_Stim*rr][q]/100 // RGB LEDs are 0-100
    Stimulus[nP_stim*rr,nP_stim*(rr+1)-1][nLEDs-1]=NoiseArray4LEDs[p-nP_Stim*rr][q]/200 // UV LED uis 0-200 so extra div 2
endfor


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Bring stimulus to 500 Hz timebase             /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// generating output arrays
make /o/n=(nP_data_upsampled,4) Stim_upsampled = 0
setscale /p x,0,LineDuration,"s" Stim_upsampled
// upsampling stimulus array 
print "upsampling Stimulus..."
for (ll=0;ll<nLEDs;ll+=1)
    for (tt=0;tt<nTriggers-1;tt+=1)
        for (pp=0;pp<1/LineDuration;pp+=1)
            variable absolutetime = (Triggertimes[tt])*(1/LineDuration)+pp // number of 2ms steps into the stimulus
            variable relativetime = (Triggertimes[tt]-Triggertimes[0])/LineDuration+pp // number of 2ms steps into the stimulus
            variable stimposition = floor(relativetime/(timebase_s_stim/LineDuration))
            do
                if (stimposition>=nP_stim)
                    stimposition-=nP_stim
                else
                    break
                endif
            while(1)
            Stim_upsampled[absolutetime][ll]=Stimulus[stimposition][ll]          
        endfor
    endfor
endfor


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// NOW FIND EVENTS ETC /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
variable nP_kernel = nSeconds_kernel/LineDuration 
variable nP_kernel_prezero = nSeconds_kernel_prezero/LineDuration 
variable nP_kernel_baseline = nSeconds_kernel_baseline/LineDuration
variable nP_kernel_eventline = nSeconds_kernel_eventline/LineDuration  
make /o/n=(nP_kernel,nLEDs,nROIs) all_kernels = 0
make /o/n=(nLEDs,nROIs) all_kernels_SD = 0
setscale /p x,-nSeconds_kernel_prezero,LineDuration,"s" all_kernels
// find events in data
print "calculating kernels..."

for (rr=0;rr<nROIs;rr+=1)
    print "ROI", rr, "/",nRois-1
    make /o/n=(nP_data) currentwave = InputTraces[p][rr]
    smooth /DIM=0 1, currentwave  // smooth before caliculating DIF 07012017 TY
    Differentiate/DIM=0  currentwave/D=currentwave_DIF
    Wavestats/Q currentwave_DIF
    currentwave_DIF/=V_SDev // normalise to SDs
    
    for (pp=floor(Triggertimes[0]/timebase_s_data);pp<Triggertimes[nTriggers-1]/timebase_s_data;pp+=1)
		if (currentwave_DIF[pp]>Noise_Threshold)        
           		 Multithread all_kernels[][][rr]+=Stim_upsampled[(pp+1)*(timebase_s_data/LineDuration)+InputTraceTimes[0][rr]/LineDuration-nP_kernel_prezero+p][q] * currentwave_DIF[pp]  //add 1 to pp to counter shift in DIF  07012017 TY
	       endif
    endfor
    // normalise each kernel & check quality
    for (ll=0;ll<nLEDs;ll+=1)
        make /o/n=(nP_kernel_baseline) currentkernel = all_kernels[p][ll][rr]
        Wavestats/Q currentkernel
        all_kernels[][ll][rr]-=V_Avg
        all_kernels[][ll][rr]/=V_SDev
        
        make /o/n=(nP_kernel_eventline) currentkernel = all_kernels[p+nP_kernel-nP_kernel_eventline][ll][rr]
        Wavestats/Q currentkernel
        all_kernels_SD[ll][rr]=V_SDev
    endfor
    
endfor

// export handling
duplicate /o all_kernels $output_name1
// display function
if (display_kernels==1)
    display /k=1
    
    make /o/n=(4,3) RGBU_Colours = 0
    RGBU_Colours[0][0]=65535 // Red
    RGBU_Colours[1][1]=65535 // Green
    RGBU_Colours[2][2]=65535 // Blue
    RGBU_Colours[3][0]=65535/2 // UV
    RGBU_Colours[3][2]=65535/2 // UV
    
    for (rr=0;rr<nRois;rr+=1)
        string YAxisName = "YAxis_Roi"+Num2Str(rr)
        string tracename
        for (ll=0;ll<nLEDs;ll+=1)
            tracename = output_name1+"#"+Num2Str(rr*nLEDs+ll)
            if (ll==0 && rr==0)
                tracename = output_name1
            endif
            Appendtograph /l=$YAxisName $output_name1[][ll][rr]
            
            ModifyGraph rgb($tracename)=(RGBU_Colours[ll][0],RGBU_Colours[ll][1],RGBU_Colours[ll][2])
            
            if (all_kernels_SD[ll][rr]>highlightSD)
                ModifyGraph lsize($tracename)=1.5
            elseif (all_kernels_SD[ll][rr]<suppressSD)
                ModifyGraph lsize($tracename)=0.5
            endif
            
        endfor  
        
        variable plotfrom = 1-((rr+1)/nRois)
        variable plotto = 1-(rr/nRois)
        
        ModifyGraph fSize($YAxisName)=8,axisEnab($YAxisName)={plotfrom,plotto};DelayUpdate
        ModifyGraph freePos($YAxisName)={0,kwFraction};DelayUpdate
        Label $YAxisName "\\Z10"+Num2Str(rr)
        ModifyGraph noLabel($YAxisName)=1,axThick($YAxisName)=0;DelayUpdate
        ModifyGraph lblRot($YAxisName)=-90
        
       SetAxis $YAxisName -nSDplot,nSDplot
    endfor
    ModifyGraph fSize(bottom)=8,axisEnab(bottom)={0.05,1};DelayUpdate
    Label bottom "\\Z10Time (\U)"
    ModifyGraph zero(bottom)=3
endif
    
    
    
// cleanup
killwaves InputTraces, InputTraceTimes, currentkernel, currentwave,currentwave_DIF ,all_kernels, all_kernels_SD 
killwaves stim_upsampled // comment to check noise stim speed if in doubt

end

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_IPLKernels()
// 1 // check for Parameter Table
if (waveexists($"NoiseArray4LEDs")==0)
    print "Warning: NoiseArray4LEDs wave missing - please import! Procedure aborted."
    abort
endif
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
// 5 //  check if Kernels are there
if (waveexists($"kernels"+Num2Str(Channel))==0)
    print "Warning: Kernels not yet generated - doing that now..."
    OS_LED_Noise()
    DoUpdate
endif

// flags from "OS_Parameters"
variable Display_IPL = OS_Parameters[%Display_Stuff]
variable LineDuration = OS_Parameters[%LineDuration]
variable nSeconds_kernel = OS_Parameters[%Noise_FilterLength_s] // nSD over baseline in time differential //read from OS_Parameters
variable YRange = OS_Parameters[%Kernel_SDHistPlot]
variable Kernel_PixelWiseMap = OS_Parameters[%Kernel_MapPxbyPx] // 0 or 1
variable Kernel_MapSmth_um = OS_Parameters[%Kernel_MapSmth] // in microns
variable Kernel_MapRange = OS_Parameters[%Kernel_MapRange]

// data handling
if (waveexists($"Positions")==0)
    print "Warning: -Positions- wave not generated yet from SARFIA - please do that first! - aborted"
    abort
endif
wave wParamsNum // Reads data-header
wave positions
wave CoM//GeoC
wave ROIs
string kernels_name = "kernels"+Num2Str(Channel)
duplicate /o $kernels_name InputKernels
string output_name1 = "KernelMap"+Num2Str(Channel)
variable pp,ll,rr

variable nSeconds_kernel_prezero = nSeconds_kernel-0.3
variable nSeconds_kernel_baseline = 0.2
variable nSeconds_kernel_eventline = 0.8 // last X s

variable Peak1Start = 375-15 // could be more elegant
variable Peak1End = 375+15
variable Peak2Start = 464-15
variable Peak2End = 464+15
variable KernelStart = 200
variable KernelEnd = 550

/// 
// calculating basic parameters
variable nP = DimSize(InputKernels,0)
variable nLEDs = DimSize(InputKernels,1)
variable nRois = DimSize(InputKernels,2)
variable nX = Dimsize(Rois,0)
variable nY = Dimsize(Rois,1)
variable zoom = wParamsNum(30) // extract zoom
variable px_Size = (0.65/zoom * 110)/nX // microns
variable Kernel_MapSmth_px = Kernel_MapSmth_um/px_Size

variable kernel_eventline_start = nP - nSeconds_kernel_eventline / LineDuration - 1
variable kernel_eventline_nP = nSeconds_kernel_eventline / LineDuration

// evaluate kernels

make /o/n=(nX,nY,nP) KernelMap_SD_StackR = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackG = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackB = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackU = 0
make /o/n=(nRois,nLEDs) KernelList_SD = NaN
make /o/n=(nX,nY,nLEDs) KernelMap_SD = 0
make /o/n=(nRois,nLEDs) KernelList_biphasic = NaN
make /o/n=(nX,nY,nLEDs) KernelMap_biphasic = 0
make /o/n=(nRois,nLEDs) KernelList_FFT = NaN
make /o/n=(nX,nY,nLEDs) KernelMap_FFT = 0

setscale /p x,-nX/2 * px_Size, px_Size,"microns" KernelMap_SD_StackR, KernelMap_SD_StackG, KernelMap_SD_StackB, KernelMap_SD_StackU
setscale /p y,-nY/2 * px_Size, px_Size,"microns" KernelMap_SD_StackR, KernelMap_SD_StackG, KernelMap_SD_StackB, KernelMap_SD_StackU
setscale /p x,-nX/2 * px_Size, px_Size,"microns" KernelMap_SD, KernelMap_biphasic, kernelMap_FFT
setscale /p y,-nY/2 * px_Size, px_Size,"microns" KernelMap_SD, KernelMap_biphasic, kernelMap_FFT

make /o/n=(Peak1end-Peak1Start) currentwavePeak1 = 0
make /o/n=(Peak2end-Peak2Start) currentwavePeak2 = 0
make /o/n=(KernelEnd-KernelStart) currentwaveplus = 0
make /o/n=(KernelEnd-KernelStart) currentwaveminus = 0
make /o/n=(KernelEnd-KernelStart) currentwave = 0

for (rr=0;rr<nRois;rr+=1)
	for (ll=0;ll<nLEDs;ll+=1)

		 // get kernel peak SD 
		 Multithread currentwavePeak1[]=InputKernels[p+Peak1Start][ll][rr]
		 Wavestats/Q currentwavePeak1
		 variable Peak1 = V_Avg
		 Multithread currentwavePeak2[]=InputKernels[p+Peak2Start][ll][rr]	 
		 Wavestats/Q currentwavePeak2
		 variable Peak2 = V_Avg		 
  		 KernelList_SD[rr][ll]=Peak2 - Peak1
		 KernelMap_SD[CoM[rr][0]][CoM[rr][1]][ll]=Peak2 - Peak1
	
		// get kernel biphasicity index
		 Multithread currentwaveplus[]=InputKernels[p+KernelStart][ll][rr]
		 Multithread currentwaveminus[]=InputKernels[p+KernelStart][ll][rr]	
		 currentwaveplus[]=(currentwaveplus[p]<0)?(0):(currentwaveplus[p])
		 currentwaveminus[]=(currentwaveminus[p]>0)?(0):(currentwaveminus[p])		 
		 Wavestats/Q Currentwaveplus
		 variable AUCPlus = V_Sum
		 Wavestats/Q Currentwaveminus
		 variable AUCMinus = -V_Sum
		 KernelList_biphasic[rr][ll]= 1-(((AUCPlus - AUCMinus) / (AUCPlus + AUCMinus))^2)^0.5
		 KernelMap_biphasic[CoM[rr][0]][CoM[rr][1]][ll] = 1-(((AUCPlus - AUCMinus) / (AUCPlus + AUCMinus))^2)^0.5

		// get Center of Mass of the Fourier Transform (i.e. the spectral centroid)
		Multithread currentwave[]=InputKernels[p+KernelStart][ll][rr]
		setscale/p x,0,LineDuration,"s" currentwave
		FFT/OUT=4/PAD={1024}/DEST=currentwave_FFT currentwave
		Wavestats/Q currentwave_FFT
		variable ff
		variable SumOfMagTimesFrequency = 0
		variable frequencySpacing = 0.48828 // Hz - from FFT at given kernel length		
		for (ff=0;ff<20;ff+=1) // only doing the 1st 20 points (up to about 9 Hz) as there is basically no power after that
			SumOfMagTimesFrequency+=FrequencySpacing*ff*currentwave_FFT[ff]
		endfor
		variable SumOfMagnitudes = V_Sum
		variable SpectralCentroid = SumOfMagTimesFrequency/SumOfMagnitudes
		 KernelList_FFT[rr][ll]= SpectralCentroid
		 KernelMap_FFT[CoM[rr][0]][CoM[rr][1]][ll] = SpectralCentroid
	endfor
	
	if (Kernel_PixelWiseMap==1)
	 	Multithread KernelMap_SD_StackR[CoM[rr][0]][CoM[rr][1]][]=InputKernels[r][ll][rr]
	 	Multithread KernelMap_SD_StackG[CoM[rr][0]][CoM[rr][1]][]=InputKernels[r][ll][rr]
		Multithread KernelMap_SD_StackB[CoM[rr][0]][CoM[rr][1]][]=InputKernels[r][ll][rr] 	 	
		Multithread KernelMap_SD_StackU[CoM[rr][0]][CoM[rr][1]][]=InputKernels[r][ll][rr] 	 		 	
	else
	 	Multithread KernelMap_SD_StackR[][][]=(ROIs[p][q]==-rr-1)?(InputKernels[r][0][rr]):(KernelMap_SD_StackR[p][q][r])
	 	Multithread KernelMap_SD_StackG[][][]=(ROIs[p][q]==-rr-1)?(InputKernels[r][1][rr]):(KernelMap_SD_StackG[p][q][r])
 	 	Multithread KernelMap_SD_StackB[][][]=(ROIs[p][q]==-rr-1)?(InputKernels[r][2][rr]):(KernelMap_SD_StackB[p][q][r])
	 	Multithread KernelMap_SD_StackU[][][]=(ROIs[p][q]==-rr-1)?(InputKernels[r][3][rr]):(KernelMap_SD_StackU[p][q][r])
	 	Multithread KernelMap_SD[][][]=(ROIs[p][q]==-rr-1)?(KernelList_SD[rr][r]):(KernelMap_SD[p][q][r])
	 	Multithread KernelMap_biphasic[][][]=(ROIs[p][q]==-rr-1)?(KernelList_biphasic[rr][r]):(KernelMap_biphasic[p][q][r])
	 	Multithread KernelMap_FFT[][][]=(ROIs[p][q]==-rr-1)?(KernelList_FFT[rr][r]):(KernelMap_FFT[p][q][r])
	 endif
endfor
// smooth maps? (default = 0)
if (Kernel_MapSmth_px>0)
 	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_SD_StackR
 	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_SD_StackG
 	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_SD_StackB
 	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_SD_StackU	 		 		 	
 	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_SD_StackR
 	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_SD_StackG
 	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_SD_StackB
 	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_SD_StackU

	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_SD
	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_SD
	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_biphasic
	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_biphasic	
	Smooth /DIM=0 Kernel_MapSmth_px, KernelMap_FFT			
	Smooth /DIM=1 Kernel_MapSmth_px, KernelMap_FFT			
endif	 		 		 	

// Make RGB overview maps
// uses Kernel_MapRange to equalise images to 16 bits
make /o/n=(nX,nY,3) KernelMap_SD_RB = 0
KernelMap_SD_RB[][][0]=((KernelMap_SD[p][q][0]/Kernel_MapRange)+1)*2^15
//KernelMap_SD_RB[][][1]=((KernelMap_SD[p][q][0]/Kernel_MapRange)+1)*2^15
KernelMap_SD_RB[][][2]=((KernelMap_SD[p][q][2]/Kernel_MapRange)+1)*2^15
KernelMap_SD_RB[][][]=(KernelMap_SD_RB[p][q][r]>2^16-1)?(2^16-1):(KernelMap_SD_RB[p][q][r])
KernelMap_SD_RB[][][]=(KernelMap_SD_RB[p][q][r]<0)?(0):(KernelMap_SD_RB[p][q][r])
KernelMap_SD_RB[][][]=(KernelMap_SD_RB[p][q][0]==2^15 && KernelMap_SD_RB[p][q][2]==2^15)?(2^15):(KernelMap_SD_RB[p][q][r])

make /o/n=(nX,nY,3) KernelMap_SD_GU = 0
KernelMap_SD_GU[][][1]=((KernelMap_SD[p][q][1]/Kernel_MapRange)+1)*2^15
KernelMap_SD_GU[][][0]=((KernelMap_SD[p][q][3]/Kernel_MapRange)+1)*2^15
KernelMap_SD_GU[][][2]=((KernelMap_SD[p][q][3]/Kernel_MapRange)+1)*2^15
KernelMap_SD_GU[][][]=(KernelMap_SD_GU[p][q][r]>2^16-1)?(2^16-1):(KernelMap_SD_GU[p][q][r])
KernelMap_SD_GU[][][]=(KernelMap_SD_GU[p][q][r]<0)?(0):(KernelMap_SD_GU[p][q][r])

make /o/n=(nX,nY,3) KernelMap_SD_RGB = 0
KernelMap_SD_RGB[][][0]=((KernelMap_SD[p][q][0]/Kernel_MapRange)+1)*2^15
KernelMap_SD_RGB[][][1]=((KernelMap_SD[p][q][1]/Kernel_MapRange)+1)*2^15
KernelMap_SD_RGB[][][2]=((KernelMap_SD[p][q][2]/Kernel_MapRange)+1)*2^15
KernelMap_SD_RGB[][][]=(KernelMap_SD_RGB[p][q][r]>2^16-1)?(2^16-1):(KernelMap_SD_RGB[p][q][r])
KernelMap_SD_RGB[][][]=(KernelMap_SD_RGB[p][q][r]<0)?(0):(KernelMap_SD_RGB[p][q][r])

make /o/n=(nX,nY,3) KernelMap_SD_RGBU = 0
KernelMap_SD_RGBU[][][0]=((KernelMap_SD[p][q][0]/Kernel_MapRange)+1)*2^14
KernelMap_SD_RGBU[][][1]=((KernelMap_SD[p][q][1]/Kernel_MapRange)+1)*2^15
KernelMap_SD_RGBU[][][2]=((KernelMap_SD[p][q][2]/Kernel_MapRange)+1)*2^14
KernelMap_SD_RGBU[][][0]+=((KernelMap_SD[p][q][3]/Kernel_MapRange)+1)*2^14
KernelMap_SD_RGBU[][][2]+=((KernelMap_SD[p][q][3]/Kernel_MapRange)+1)*2^14
KernelMap_SD_RGBU[][][]=(KernelMap_SD_RGBU[p][q][r]>2^16-1)?(2^16-1):(KernelMap_SD_RGBU[p][q][r])
KernelMap_SD_RGBU[][][]=(KernelMap_SD_RGBU[p][q][r]<0)?(0):(KernelMap_SD_RGBU[p][q][r])

setscale /p x,-nX/2 * px_Size, px_Size,"microns" KernelMap_SD_RGBU, KernelMap_SD_RGB, KernelMap_SD_RB, KernelMap_SD_GU
setscale /p y,-nY/2 * px_Size, px_Size,"microns" KernelMap_SD_RGBU, KernelMap_SD_RGB, KernelMap_SD_RB, KernelMap_SD_GU
	
// get IPL depth profiles
variable nBins = 25

make /o/n=(nBins,nLEDs) kernel_IPLHists = 0
make /o/n=(nBins) nPixels_per_depth = 0

for (rr=0;rr<nRois;rr+=1)
	for (ll=0;ll<nLEDs;ll+=1)
		kernel_IPLHists[floor(positions[rr]/(100/nBins))][ll]+=KernelList_SD[rr][ll]
		nPixels_per_depth[floor(positions[rr]/(100/nBins))]+=1
	endfor
endfor


setscale x,0,100,"%" kernel_IPLHists

// display
if (display_IPL == 1)
	Display /k=1
	Appendimage /l=Y1 /b=X1 KernelMap_SD_RB
	Appendimage /l=Y1 /b=X2 KernelMap_SD_GU
	Appendimage /l=Y2 /b=X1 KernelMap_SD_RGB
	Appendimage /l=Y2 /b=X2 KernelMap_SD_RGBU
	ModifyGraph fSize=8,lblPos=47,axisEnab(Y1)={0.51,1},axisEnab(X1)={0,0.49};DelayUpdate
	ModifyGraph axisEnab(X2)={0.51,1},axisEnab(Y2)={0,0.49};DelayUpdate
	ModifyGraph freePos(Y1)={0,kwFraction},freePos(X1)={0,kwFraction};DelayUpdate
	ModifyGraph freePos(X2)={0,kwFraction},freePos(Y2)={0,kwFraction}
	
	
//	display /k=1
//	Appendtograph /l=RedY kernel_IPLHists[][0]
//	Appendtograph /l=GreenY kernel_IPLHists[][1]
//	Appendtograph /l=BlueY kernel_IPLHists[][2]
//	Appendtograph /l=UVY kernel_IPLHists[][3]
//	
//	ModifyGraph fSize=8,axisEnab(bottom)={0.05,1};DelayUpdate
//	Label bottom "\\Z10IPL depth (\U)";DelayUpdate
//	SetAxis bottom 0,100
//
//	ModifyGraph axisEnab(UVY)={0.05,0.25},axisEnab(BlueY)={0.3,0.5};DelayUpdate
//	ModifyGraph axisEnab(GreenY)={0.55,0.75},axisEnab(RedY)={0.8,1};DelayUpdate
//	ModifyGraph freePos(RedY)={0,kwFraction},freePos(GreenY)={0,kwFraction};DelayUpdate
//	ModifyGraph freePos(BlueY)={0,kwFraction},freePos(UVY)={0,kwFraction};DelayUpdate
//	SetAxis RedY -YRange,YRange;DelayUpdate
//	SetAxis GreenY -YRange,YRange;DelayUpdate
//	SetAxis BlueY -YRange,YRange;DelayUpdate
//	SetAxis UVY -YRange,YRange
//
//
//	ModifyGraph axThick(RedY)=1,lblPos(RedY)=47, nticks(RedY)=2
//	ModifyGraph axThick(GreenY)=1,lblPos(GreenY)=47, nticks(GreenY)=2
//	ModifyGraph axThick(BlueY)=1,lblPos(BlueY)=47, nticks(BlueY)=2
//	ModifyGraph axThick(UVY)=1,lblPos(UVY)=47, nticks(UVY)=2			
//	Label RedY "\\Z10Peak-to-Peak SD (R)"
//	Label GreenY "\\Z10Peak-to-Peak SD (G)"
//	Label BlueY "\\Z10Peak-to-Peak SD (B)"
//	Label UVY "\\Z10Peak-to-Peak SD (U)"			
//
//	ModifyGraph rgb(kernel_IPLHists#1)=(0,52224,0)
//	ModifyGraph rgb(kernel_IPLHists#2)=(0,0,65280)
//	ModifyGraph rgb(kernel_IPLHists#3)=(29440,0,58880)
//	ModifyGraph mode=5, lsize = 1.5, hbFill=5
//	
//	•Display /k=1 kernel_IPLHists[][0],kernel_IPLHists[][1],kernel_IPLHists[][2],kernel_IPLHists[][3]
//•ModifyGraph mode=6,lsize=1.5,rgb(kernel_IPLHists#1)=(0,52224,0);DelayUpdate
//•ModifyGraph rgb(kernel_IPLHists#2)=(0,0,65280);DelayUpdate
//•ModifyGraph rgb(kernel_IPLHists#3)=(29440,0,58880)
//•ModifyGraph zero(left)=2,fSize=8,axisEnab(left)={0.05,1};DelayUpdate
//•ModifyGraph axisEnab(bottom)={0.05,1};DelayUpdate
//•Label left "\\Z10Peak-to-Peak SD";DelayUpdate
//•Label bottom "\\Z10IPL depth (\\U)"
//•ModifyGraph mode=8,toMode=1
//ModifyGraph marker=19,msize=3
	
	
endif


// cleanup
killwaves InputKernels,nPixels_per_depth, currentwavePeak1, currentwavePeak2, currentwaveplus, currentwaveminus, currentwave,currentwave_FFT

end