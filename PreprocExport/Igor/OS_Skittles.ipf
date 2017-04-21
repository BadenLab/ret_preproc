#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_Skittles()

// 1 // check for Parameter Table
if (waveexists($"SkittlesNoise")==0)
    print "Warning: SkittlesNoise wave missing - please import! Procedure aborted."
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
wave SkittlesNoise // official stimulus array
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

// calculating basic parameters
variable nP_stim = Dimsize(SkittlesNoise,0)
variable nP_data = Dimsize(InputTraceTimes,0)
variable nLEDs = Dimsize(SkittlesNoise,1)
variable nTriggers = Dimsize(Triggertimes,0)
variable timebase_s_stim =noise_interval
variable timebase_s_data = InputTraceTimes[1][0]-InputTraceTimes[0][0]
variable nP_data_upsampled = ceil(nP_data * timebase_s_data * 1/LineDuration)
variable nP_stim_upsampled = ceil(nP_stim * timebase_s_stim * 1/LineDuration)
variable nStim_repeats = ceil(nP_data_upsampled / nP_stim_upsampled )
make /o/n=(nP_stim*nStim_repeats,nLEDs) Stimulus = NaN
for (rr=0;rr<nStim_repeats;rr+=1)
    Stimulus[nP_stim*rr,nP_stim*(rr+1)-1][0,nLEDs-1]=SkittlesNoise[p-nP_Stim*rr][q] 
endfor


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Bring stimulus to 500 Hz timebase             /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// generating output arrays
make /o/n=(nP_data_upsampled,nLEDs) Stim_upsampled = 0
setscale /p x,0,LineDuration,"s" Stim_upsampled
// upsampling stimulus array 
print "upsampling Stimulus..."

variable LoopDuration = timebase_s_stim * nP_stim


for (tt=1;tt<nTriggers;tt+=1) // note starting from Trigger 1 not 0
	variable absolutetime = (Triggertimes[tt])*(1/LineDuration)+pp - (timebase_s_stim/LineDuration)/2 // number of 2ms steps into the stimulus 
	Stim_upsampled[absolutetime,absolutetime+LoopDuration/LineDuration][]=Stimulus[(p-absolutetime)*LineDuration/timebase_s_stim][q]
       
endfor



///////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// NOW FIND EVENTS ETC /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    
    make /o/n=(nLEDs,3) RGBU_Colours = 0
    make /o/n=(1) M_colors
    ColorTab2Wave Rainbow256
    for (ll=0;ll<nLEDs;ll+=1)
	 	RGBU_Colours[ll][]=M_Colors[256/nLEDs * ll][q]
	
	
    endfor


    
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
//killwaves stim_upsampled // comment to check noise stim speed if in doubt

print "to display individual kernels, call OS_PlotKernels(Roinumber)"

end


end