#pragma rtGlobals=3		// Use modern global access method and strict wave access.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// SKITTLES SWEEP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_SkittlesSweep()


variable nLEDs = 19
make /o/n=(nLEDs) SkittlesWavelengths = {671,641,615,598,572,557,535,519,505,494,480,466,446,424,407,393,368,356,320}
variable ReadoutTimes_s = 0.2 // i.e. 100 ms after Start and before End of step
variable ReadoutWindow_s = 0.3 // i.e. integrate 100 ms worth of trace


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
// 5 //  check if Averages"N" is there
if (waveexists($"Averages"+Num2Str(Channel))==0)
	print "Warning: Averages wave not yet generated - doing that now..."
	OS_BasicAveraging()
	DoUpdate
endif

// flags from "OS_Parameters"
variable LineDuration = OS_Parameters[%LineDuration]
variable Display_Tunings = OS_Parameters[%Display_Stuff]

// data handling
string input_name = "Snippets"+Num2Str(Channel)
duplicate /o $input_name InputData

variable nP = DimSize(InputData,0)
variable nFlashes = DimSize(InputData,1)
variable nRois = DimSize(InputData,2)

string output_name1 = "SweepMapMeans"+Num2Str(Channel)
string output_name2 = "SweepMapSnippets"+Num2Str(Channel)
string output_name3 = "SweepTuningMeans"+Num2Str(Channel)
string output_name4 = "SweepTuningSnippets"+Num2Str(Channel)

variable rr,ll,ff,cc

// make new Looped Snippets arrays
variable nCompleteLoops = floor(nFlashes/nLEDs)
print nCompleteLoops, "complete loops of", nLEDs, "LEDs"

make /o/n=(nP,nLEDs,nROIs) SweepMeans=0
make /o/n=(nP,nLEDs,nCompleteLoops,nROIs) SweepSnippets=0

setscale /p x,0,LineDuration,"s" SweepMeans,SweepSnippets

variable CurrentLED = 0
variable CurrentLoop = 0
for (ff=0;ff<nFlashes;ff+=1)
	SweepSnippets[][CurrentLED][CurrentLoop][]=InputData[p][ff][s]
	SweepMeans[][CurrentLED][]+=InputData[p][ff][r]/nCompleteLoops
	CurrentLED+=1
	if (CurrentLED>=nLEDs)
		CurrentLED=0
		CurrentLoop+=1
	endif
endfor

//Extract Tunings

variable ONPeakTime_P = ReadoutTimes_s/LineDuration
variable ONsusTime_P = nP/2 - ReadoutTimes_s/LineDuration
variable OFFPeakTime_P = nP/2 + ReadoutTimes_s/LineDuration
variable OFFsustime_P = nP - ReadoutTimes_s/LineDuration

make /o/n=(nLEDs,4,nROIs) SweepTuning_mean = NaN
make /o/n=(nLEDs,4,nCompleteLoops,nROIs) SweepTuning_snippets = NaN
make /o/n=(ReadoutWindow_s/LineDuration) currentwave = 0

for (rr=0;rr<nROIs;rr+=1)
	for (ll=0;ll<nLEDs;ll+=1)
	
		Multithread currentwave[]=SweepMeans[p+ONPeakTime_P][ll][rr]
		Wavestats/Q currentwave
		SweepTuning_mean[ll][0][rr]=V_Avg
		
		Multithread currentwave[]=SweepMeans[-p+ONsusTime_P][ll][rr]
		Wavestats/Q currentwave
		SweepTuning_mean[ll][1][rr]=V_Avg

		Multithread currentwave[]=SweepMeans[p+OFFPeakTime_P][ll][rr]
		Wavestats/Q currentwave
		SweepTuning_mean[ll][2][rr]=V_Avg

		Multithread currentwave[]=SweepMeans[-p+OFFsusTime_P][ll][rr]
		Wavestats/Q currentwave
		SweepTuning_mean[ll][3][rr]=V_Avg

		for (cc=0;cc<nCompleteLoops;cc+=1)

			Multithread currentwave[]=SweepSnippets[p+ONPeakTime_P][ll][cc][rr]
			Wavestats/Q currentwave
			SweepTuning_Snippets[ll][0][cc][rr]=V_Avg
			
			Multithread currentwave[]=SweepSnippets[p+ONsusTime_P][ll][cc][rr]
			Wavestats/Q currentwave
			SweepTuning_Snippets[ll][1][cc][rr]=V_Avg
			
			Multithread currentwave[]=SweepSnippets[p+OFFPeakTime_P][ll][cc][rr]
			Wavestats/Q currentwave
			SweepTuning_Snippets[ll][2][cc][rr]=V_Avg
			
			Multithread currentwave[]=SweepSnippets[p+OFFsusTime_P][ll][cc][rr]
			Wavestats/Q currentwave
			SweepTuning_Snippets[ll][3][cc][rr]=V_Avg

		endfor

	endfor
endfor

// export handling
duplicate /o SweepMeans $output_name1
duplicate /o SweepSnippets $output_name2
duplicate /o SweepTuning_mean $output_name3
duplicate /o SweepTuning_snippets $output_name4

// display

if (display_tunings==1)

	display /k=1
	make /o/n=(1) M_Colors
	Colortab2Wave Rainbow256
	
	for (rr=0;rr<nRois;rr+=1)
		string YAxisName = "YAxis_Roi"+Num2Str(rr)
		string tracename
		for (ll=0;ll<nCompleteLoops;ll+=1)
			tracename = output_name4+"#"+Num2Str((rr*nCompleteLoops+ll)*4)
			if (ll==0 && rr==0)
				tracename = output_name4
			endif
			Appendtograph /l=$YAxisName /b=XOnTr $output_name4[][0][ll][rr] vs SkittlesWavelengths // ON transient
			ModifyGraph rgb($tracename)=(52224,52224,52224)
			
			tracename = output_name4+"#"+Num2Str((rr*nCompleteLoops+ll)*4+1)	
			Appendtograph /l=$YAxisName /b=XONsus $output_name4[][1][ll][rr] vs SkittlesWavelengths // ON sustained
			ModifyGraph rgb($tracename)=(52224,52224,52224)
			
			tracename = output_name4+"#"+Num2Str((rr*nCompleteLoops+ll)*4+2)	
			Appendtograph /l=$YAxisName /b=XOffTr $output_name4[][2][ll][rr] vs SkittlesWavelengths // OFF transient
			ModifyGraph rgb($tracename)=(52224,52224,52224)
			
			tracename = output_name4+"#"+Num2Str((rr*nCompleteLoops+ll)*4+3)	
			Appendtograph /l=$YAxisName /b=XOffSus $output_name4[][3][ll][rr] vs SkittlesWavelengths // OFF sustained
			ModifyGraph rgb($tracename)=(52224,52224,52224)
			
		endfor	
		
		tracename = output_name3+"#"+Num2Str(rr*4)
		if (rr==0)
			tracename = output_name3
		endif
		Appendtograph /l=$YAxisName /b=XOnTr $output_name3[][0][rr] vs SkittlesWavelengths // ON tr Means
		variable colorposition = 255 * (rr+1)/nRois
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
		ModifyGraph lsize($tracename)=1.5
		
		tracename = output_name3+"#"+Num2Str(rr*4+1)
		Appendtograph /l=$YAxisName /b=XOnSus $output_name3[][1][rr] vs SkittlesWavelengths // ON sus Means
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
		ModifyGraph lsize($tracename)=1.5

		tracename = output_name3+"#"+Num2Str(rr*4+2)		
		Appendtograph /l=$YAxisName /b=XOfftr $output_name3[][2][rr] vs SkittlesWavelengths // OFF tr Means
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
		ModifyGraph lsize($tracename)=1.5
		
		tracename = output_name3+"#"+Num2Str(rr*4+3)
		Appendtograph /l=$YAxisName /b=XOffSus $output_name3[][3][rr] vs SkittlesWavelengths // OFF sus Means
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
		ModifyGraph lsize($tracename)=1.5
		
		variable plotfrom = (1-((rr+1)/nRois))*0.85+0.05
		variable plotto = (1-(rr/nRois))*0.85+0.05
		
		ModifyGraph fSize($YAxisName)=8,axisEnab($YAxisName)={plotfrom,plotto};DelayUpdate
		ModifyGraph freePos($YAxisName)={0,kwFraction};DelayUpdate
		Label $YAxisName "\\Z10"+Num2Str(rr)
		ModifyGraph noLabel($YAxisName)=1,axThick($YAxisName)=0;DelayUpdate
		ModifyGraph lblRot($YAxisName)=-90
	endfor
	
	ModifyGraph fSize=8,lblPos(XOnTr)=47,axisEnab(XOnTr)={0.05,0.25};DelayUpdate
	ModifyGraph freePos(XOnTr)={0,kwFraction}
	Label XOnTr "\\Z10Wavelength (nm)"
	SetAxis XOnTr 300,700

	ModifyGraph fSize=8,lblPos(XOnSus)=47,axisEnab(XOnSus)={0.3,0.5};DelayUpdate
	ModifyGraph freePos(XOnSus)={0,kwFraction}
	Label XOnSus "\\Z10Wavelength (nm)"
	SetAxis XOnSus 300,700
	
	ModifyGraph fSize=8,lblPos(XOffTr)=47,axisEnab(XOffTr)={0.55,0.75};DelayUpdate
	ModifyGraph freePos(XOffTr)={0,kwFraction}
	Label XOffTr "\\Z10Wavelength (nm)"
	SetAxis XOffTr 300,700
	
	ModifyGraph fSize=8,lblPos(XOffSus)=47,axisEnab(XOffSus)={0.8,1};DelayUpdate
	ModifyGraph freePos(XOffSus)={0,kwFraction}
	Label XOffSus "\\Z10Wavelength (nm)"
	SetAxis XOffSus 300,700

	•ShowTools/A arrow
	•SetDrawEnv xcoord= XOnTr,fstyle= 1, fsize= 10;DelayUpdate
	•DrawText 360,0.025,"ON transient"
	•SetDrawEnv xcoord= XOffTr,fstyle= 1,  fsize= 10;DelayUpdate
	•DrawText 360,0.025,"OFF transient"
	•SetDrawEnv xcoord= XOnSus,fstyle= 1, fsize= 10;DelayUpdate
	•DrawText 360,0.025,"ON sustained"
	•SetDrawEnv xcoord= XOffSus,fstyle= 1,  fsize= 10;DelayUpdate
	•DrawText 360,0.025,"OFF sustained"
	HideTools/A
	

endif


// cleanup
killwaves InputData, currentwave


end

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// OPSIN PLOTTER
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_SkittlesSweep_Plot(rr)
variable rr

wave SweepTuningMeans0
wave SweepTuningSnippets0
wave SkittlesWavelengths
wave zf_opsins
variable ll
variable nCompleteLoops = Dimsize(SweepTuningSnippets0,2)

display /k=1
	
string YAxisName = "YAxis"
string tracename
	for (ll=0;ll<nCompleteLoops;ll+=1)
		tracename = "SweepTuningSnippets0#"+Num2Str((ll)*4)
		if (ll==0 && rr==0)
			tracename = "SweepTuningSnippets0"
		endif
		Appendtograph /l=$YAxisName /b=XOnTr SweepTuningSnippets0[][0][ll][rr] vs SkittlesWavelengths // ON transient
		ModifyGraph rgb($tracename)=(52224,52224,52224)
		
		tracename = "SweepTuningSnippets0#"+Num2Str((ll)*4+1)	
		Appendtograph /l=$YAxisName /b=XONsus SweepTuningSnippets0[][1][ll][rr] vs SkittlesWavelengths // ON sustained
		ModifyGraph rgb($tracename)=(52224,52224,52224)
		
		tracename = "SweepTuningSnippets0#"+Num2Str((ll)*4+2)	
		Appendtograph /l=$YAxisName /b=XOffTr SweepTuningSnippets0[][2][ll][rr] vs SkittlesWavelengths // OFF transient
		ModifyGraph rgb($tracename)=(52224,52224,52224)
		
		tracename = "SweepTuningSnippets0#"+Num2Str((ll)*4+3)	
		Appendtograph /l=$YAxisName /b=XOffSus SweepTuningSnippets0[][3][ll][rr] vs SkittlesWavelengths // OFF sustained
		ModifyGraph rgb($tracename)=(52224,52224,52224)
		
	endfor	
		

	tracename = "SweepTuningMeans0"
	Appendtograph /l=$YAxisName /b=XOnTr SweepTuningMeans0[][0][rr] vs SkittlesWavelengths // ON tr Means
	ModifyGraph rgb($tracename)=(0,0,0)
	ModifyGraph lsize($tracename)=1.5
	
	tracename = "SweepTuningMeans0#"+Num2Str(1)
	Appendtograph /l=$YAxisName /b=XOnSus SweepTuningMeans0[][1][rr] vs SkittlesWavelengths // ON sus Means
	ModifyGraph rgb($tracename)=(0,0,0)
	ModifyGraph lsize($tracename)=1.5

	tracename = "SweepTuningMeans0#"+Num2Str(2)		
	Appendtograph /l=$YAxisName /b=XOfftr SweepTuningMeans0[][2][rr] vs SkittlesWavelengths // OFF tr Means
	ModifyGraph rgb($tracename)=(0,0,0)
	ModifyGraph lsize($tracename)=1.5
	
	tracename = "SweepTuningMeans0#"+Num2Str(3)
	Appendtograph /l=$YAxisName /b=XOffSus SweepTuningMeans0[][3][rr] vs SkittlesWavelengths // OFF sus Means
	ModifyGraph rgb($tracename)=(0,0,0)
	ModifyGraph lsize($tracename)=1.5
		
	variable plotfrom = 0.15
	variable plotto = 1
		
	ModifyGraph fSize($YAxisName)=8,axisEnab($YAxisName)={plotfrom,plotto};DelayUpdate
	ModifyGraph freePos($YAxisName)={0,kwFraction};DelayUpdate
	Label $YAxisName "\\Z10"+Num2Str(rr)
	ModifyGraph noLabel($YAxisName)=1,axThick($YAxisName)=0;DelayUpdate
	ModifyGraph lblRot($YAxisName)=-90

///
	
	ModifyGraph fSize=8,lblPos(XOnTr)=47,axisEnab(XOnTr)={0.05,0.25};DelayUpdate
	ModifyGraph freePos(XOnTr)={0,kwFraction}
	Label XOnTr "\\Z10Wavelength (nm)"
	SetAxis XOnTr 300,700

	ModifyGraph fSize=8,lblPos(XOnSus)=47,axisEnab(XOnSus)={0.3,0.5};DelayUpdate
	ModifyGraph freePos(XOnSus)={0,kwFraction}
	Label XOnSus "\\Z10Wavelength (nm)"
	SetAxis XOnSus 300,700
	
	ModifyGraph fSize=8,lblPos(XOffTr)=47,axisEnab(XOffTr)={0.55,0.75};DelayUpdate
	ModifyGraph freePos(XOffTr)={0,kwFraction}
	Label XOffTr "\\Z10Wavelength (nm)"
	SetAxis XOffTr 300,700
	
	ModifyGraph fSize=8,lblPos(XOffSus)=47,axisEnab(XOffSus)={0.8,1};DelayUpdate
	ModifyGraph freePos(XOffSus)={0,kwFraction}
	Label XOffSus "\\Z10Wavelength (nm)"
	SetAxis XOffSus 300,700

	•ShowTools/A arrow
	•SetDrawEnv xcoord= XOnTr,fstyle= 1, fsize= 10;DelayUpdate
	•DrawText 360,0.025,"ON transient"
	•SetDrawEnv xcoord= XOffTr,fstyle= 1,  fsize= 10;DelayUpdate
	•DrawText 360,0.025,"OFF transient"
	•SetDrawEnv xcoord= XOnSus,fstyle= 1, fsize= 10;DelayUpdate
	•DrawText 360,0.025,"ON sustained"
	•SetDrawEnv xcoord= XOffSus,fstyle= 1,  fsize= 10;DelayUpdate
	•DrawText 360,0.025,"OFF sustained"
	HideTools/A
	
	//
	
	•Appendtograph /l=OpsinsY /b=XOnTr zf_opsins[][0],zf_opsins[][1],zf_opsins[][2],zf_opsins[][3]
	•Appendtograph /l=OpsinsY /b=XOnSus zf_opsins[][0],zf_opsins[][1],zf_opsins[][2],zf_opsins[][3]
	•Appendtograph /l=OpsinsY /b=XOffTr zf_opsins[][0],zf_opsins[][1],zf_opsins[][2],zf_opsins[][3]
	•Appendtograph /l=OpsinsY /b=XOffSus zf_opsins[][0],zf_opsins[][1],zf_opsins[][2],zf_opsins[][3]	
	
	•ModifyGraph fSize=8,noLabel(OpsinsY)=2,axThick(OpsinsY)=0;DelayUpdate
	•ModifyGraph axisEnab(OpsinsY)={0.05,0.3},freePos(OpsinsY)={0,kwFraction}
	
	•ModifyGraph rgb(zf_opsins)=(29440,0,58880),rgb(zf_opsins#1)=(0,0,65280);DelayUpdate
•ModifyGraph rgb(zf_opsins#2)=(0,52224,26368),rgb(zf_opsins#4)=(29440,0,58880);DelayUpdate
•ModifyGraph rgb(zf_opsins#5)=(0,0,65280),rgb(zf_opsins#6)=(0,52224,26368);DelayUpdate
•ModifyGraph rgb(zf_opsins#8)=(29440,0,58880),rgb(zf_opsins#9)=(0,0,65280);DelayUpdate
•ModifyGraph rgb(zf_opsins#10)=(0,52224,26368),rgb(zf_opsins#12)=(29440,0,58880);DelayUpdate
•ModifyGraph rgb(zf_opsins#13)=(0,0,65280),rgb(zf_opsins#14)=(0,52224,26368)
	
end


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// NOISE
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_SkittlesNoise()

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