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

// data handling
if (waveexists($"Positions")==0)
    print "Warning: -Positions- wave not generated yet from SARFIA - please do that first! - aborted"
    abort
endif
wave positions
wave GeoC
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


/// 
// calculating basic parameters
variable nP = DimSize(InputKernels,0)
variable nLEDs = DimSize(InputKernels,1)
variable nRois = DimSize(InputKernels,2)
variable nX = Dimsize(Rois,0)
variable nY = Dimsize(Rois,1)


variable kernel_eventline_start = nP - nSeconds_kernel_eventline / LineDuration - 1
variable kernel_eventline_nP = nSeconds_kernel_eventline / LineDuration

// evaluate kernels


make /o/n=(nX,nY,nLEDs*2) KernelMap_Peaks = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackR = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackG = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackB = 0
make /o/n=(nX,nY,nP) KernelMap_SD_StackU = 0
make /o/n=(nRois,nLEDs) KernelList_SD = NaN

make /o/n=(kernel_eventline_nP) currentwave = 0
make /o/n=(Peak1end-Peak1Start) currentwavePeak1 = 0
make /o/n=(Peak2end-Peak2Start) currentwavePeak2 = 0
for (rr=0;rr<nRois;rr+=1)
	for (ll=0;ll<nLEDs;ll+=1)
		 Multithread currentwave[]=InputKernels[p+kernel_eventline_start][ll][rr]
		 Wavestats/Q currentwave
		 KernelList_SD[rr][ll]=V_SDev			 
		 Multithread currentwavePeak1[]=InputKernels[p+Peak1Start][ll][rr]
		 Wavestats/Q currentwavePeak1
		 KernelMap_Peaks[GeoC[rr][0]][GeoC[rr][1]][ll*2]=V_Avg
		 Multithread currentwavePeak2[]=InputKernels[p+Peak2Start][ll][rr]	 
		 Wavestats/Q currentwavePeak2
		 KernelMap_Peaks[GeoC[rr][0]][GeoC[rr][1]][ll*2+1]=V_Avg	
	endfor
	 Multithread KernelMap_SD_StackR[GeoC[rr][0]][GeoC[rr][1]][]=InputKernels[r][0][rr]
	 Multithread KernelMap_SD_StackG[GeoC[rr][0]][GeoC[rr][1]][]=InputKernels[r][1][rr]
	 Multithread KernelMap_SD_StackB[GeoC[rr][0]][GeoC[rr][1]][]=InputKernels[r][2][rr]
	 Multithread KernelMap_SD_StackU[GeoC[rr][0]][GeoC[rr][1]][]=InputKernels[r][3][rr]
endfor


// get IPL depth profiles
variable nBins = 25

make /o/n=(nBins,nLEDs) kernel_IPLPeak1Hists = 0
make /o/n=(nBins,nLEDs) kernel_IPLPeak2Hists = 0
make /o/n=(nBins,nLEDs) kernel_IPLHists = 0
make /o/n=(nBins) nPixels_per_depth = 0

for (rr=0;rr<nRois;rr+=1)
	for (ll=0;ll<nLEDs;ll+=1)
		kernel_IPLPeak1Hists[floor(positions[rr]/(100/nBins))][ll]+=KernelMap_Peaks[GeoC[rr][0]][GeoC[rr][1]][ll*2]
		kernel_IPLPeak2Hists[floor(positions[rr]/(100/nBins))][ll]+=KernelMap_Peaks[GeoC[rr][0]][GeoC[rr][1]][ll*2+1]		
		nPixels_per_depth[floor(positions[rr]/(100/nBins))]+=1
	endfor
endfor
kernel_IPLPeak1Hists[][]/=nPixels_per_depth[p]
kernel_IPLPeak2Hists[][]/=nPixels_per_depth[p]
kernel_IPLHists[][]=kernel_IPLPeak2Hists[p][q]-kernel_IPLPeak1Hists[p][q]


setscale x,0,100,"%" kernel_IPLHists

// display
if (display_IPL == 1)
	display /k=1
	Appendtograph /l=RedY kernel_IPLHists[][0]
	Appendtograph /l=GreenY kernel_IPLHists[][1]
	Appendtograph /l=BlueY kernel_IPLHists[][2]
	Appendtograph /l=UVY kernel_IPLHists[][3]
	
	ModifyGraph fSize=8,axisEnab(bottom)={0.05,1};DelayUpdate
	Label bottom "\\Z10IPL depth (\U)";DelayUpdate
	SetAxis bottom 0,100

	ModifyGraph axisEnab(UVY)={0.05,0.25},axisEnab(BlueY)={0.3,0.5};DelayUpdate
	ModifyGraph axisEnab(GreenY)={0.55,0.75},axisEnab(RedY)={0.8,1};DelayUpdate
	ModifyGraph freePos(RedY)={0,kwFraction},freePos(GreenY)={0,kwFraction};DelayUpdate
	ModifyGraph freePos(BlueY)={0,kwFraction},freePos(UVY)={0,kwFraction};DelayUpdate
	SetAxis RedY -YRange,YRange;DelayUpdate
	SetAxis GreenY -YRange,YRange;DelayUpdate
	SetAxis BlueY -YRange,YRange;DelayUpdate
	SetAxis UVY -YRange,YRange


	ModifyGraph axThick(RedY)=1,lblPos(RedY)=47, nticks(RedY)=2
	ModifyGraph axThick(GreenY)=1,lblPos(GreenY)=47, nticks(GreenY)=2
	ModifyGraph axThick(BlueY)=1,lblPos(BlueY)=47, nticks(BlueY)=2
	ModifyGraph axThick(UVY)=1,lblPos(UVY)=47, nticks(UVY)=2			
	Label RedY "\\Z10Peak-to-Peak SD (R)"
	Label GreenY "\\Z10Peak-to-Peak SD (G)"
	Label BlueY "\\Z10Peak-to-Peak SD (B)"
	Label UVY "\\Z10Peak-to-Peak SD (U)"			

	ModifyGraph rgb(kernel_IPLHists#1)=(0,52224,0)
	ModifyGraph rgb(kernel_IPLHists#2)=(0,0,65280)
	ModifyGraph rgb(kernel_IPLHists#3)=(29440,0,58880)
	ModifyGraph mode=5, lsize = 1.5, hbFill=5
	
	•Display /k=1 kernel_IPLHists[][0],kernel_IPLHists[][1],kernel_IPLHists[][2],kernel_IPLHists[][3]
•ModifyGraph mode=6,lsize=1.5,rgb(kernel_IPLHists#1)=(0,52224,0);DelayUpdate
•ModifyGraph rgb(kernel_IPLHists#2)=(0,0,65280);DelayUpdate
•ModifyGraph rgb(kernel_IPLHists#3)=(29440,0,58880)
•ModifyGraph zero(left)=2,fSize=8,axisEnab(left)={0.05,1};DelayUpdate
•ModifyGraph axisEnab(bottom)={0.05,1};DelayUpdate
•Label left "\\Z10Peak-to-Peak SD";DelayUpdate
•Label bottom "\\Z10IPL depth (\\U)"
•ModifyGraph mode=8,toMode=1
ModifyGraph marker=19,msize=3
	
	
endif


// cleanup
killwaves InputKernels,currentwave,nPixels_per_depth

end