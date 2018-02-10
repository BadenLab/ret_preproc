#pragma rtGlobals=3		// Use modern global access method and strict wave access.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// SWOOSH MAP ////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_SkittlesSwooshMap()

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

// Hardcoded SkittleSwoosh positions - saves having to load the array
make /o/n=(80) SwooshPositionsX1 = {384,382,380,377,375,372,369,366,364,361,358,355,352,349,346,343,340,338,335,332,386,384,381,379,376,373,370,367,364,361,358,355,352,349,345,342,339,336,333,331,388,386,383,381,378,375,372,369,365,362,358,355,351,348,345,341,338,335,332,329,391,388,386,383,380,377,373,370,366,363,359,355,351,347,344,340,337,333,330,327}
make /o/n=(80) SwooshPositionsX2 = {394,391,388,385,382,379,375,371,367,363,359,355,351,347,343,339,335,331,328,325,397,394,391,388,385,381,377,373,369,364,360,355,351,346,341,337,333,329,325,322,400,398,395,391,388,384,380,376,371,366,361,355,350,345,340,335,330,326,322,318,404,402,399,395,392,388,383,378,373,368,362,356,350,344,338,332,327,323,318,314}
make /o/n=(80) SwooshPositionsY1 = {420,419,418,417,417,416,415,415,414,414,414,414,414,414,414,414,415,415,416,417,418,416,415,414,414,413,412,412,411,411,411,410,410,411,411,411,412,412,413,414,415,414,413,412,411,410,409,408,408,408,407,407,407,407,408,408,409,409,410,411,412,411,410,409,408,407,406,405,405,404,404,404,404,404,405,405,406,406,407,408}
make /o/n=(80) SwooshPositionsY2 = {410,408,407,406,405,404,403,402,402,401,401,401,401,401,402,402,403,403,404,405,407,406,405,403,402,401,400,399,399,398,398,398,398,398,398,399,400,401,402,403,405,404,402,401,399,398,397,396,396,395,395,395,395,395,395,396,397,398,399,400,403,401,400,398,397,396,394,393,393,392,392,391,391,392,392,393,394,395,396,398}
make /o/n=(160,2) SwooshPositions = NaN
SwooshPositions[0,79][1] = SwooshPositionsX1[p]
SwooshPositions[0,79][0] = SwooshPositionsY1[p]
SwooshPositions[80,159][1] = SwooshPositionsX2[p-80]
SwooshPositions[80,159][0] = SwooshPositionsY2[p-80]
killwaves SwooshPositionsX1, SwooshPositionsY1,SwooshPositionsX2, SwooshPositionsY2

// flags from "OS_Parameters"
variable Display_maps = OS_Parameters[%Display_Stuff]
variable LineDuration = OS_Parameters[%LineDuration]

// data handling
string input_name = "Averages"+Num2Str(Channel)
duplicate /o $input_name InputData

wave QualityCriterion
variable nP = DimSize(InputData,0)
variable nRois = DimSize(InputData,1)
variable nPositions = Dimsize(SwooshPositions,0)
make /o/n=(nPositions) currentwave = SwooshPositions[p][1]
WaveStats/Q currentwave
Variable MapX_min = V_Min
Variable MapX_max = V_Max
make /o/n=(nPositions) currentwave = SwooshPositions[p][0]
WaveStats/Q currentwave
Variable MapY_min = V_Min
Variable MapY_max = V_Max

string output_name1 = "SwooshMaps"+Num2Str(Channel)

variable pp,rr,cc
variable Spotdur_s = 0.2 // how long is each spot presented for (in s)
variable MapPxSize = 4
variable Smoothingfactor = 3200
variable ResponseDelay_s =  -0.0 // negative, how much do we have to bring the responses forward to compensate for any delays in visual response)
variable ResponseDelay_p = ResponseDelay_s/LineDuration

// Shift InputData arra by ResponseDelay_p points
duplicate /o inputData InputData_shift
Multithread InputData_shift[ResponseDelay_p,nP][]=InputData[p-ResponseDelay_p][q]
Multithread InputData_shift[0,ResponseDelay_p-1][]=InputData[nP-ResponseDelay_p+p][q]

// Differentiate and Smooth Input
duplicate /o InputData_shift InputData_Smth
Smooth Smoothingfactor, InputData, InputData_Smth
Differentiate/DIM=0  InputData_Smth/D=InputData_Smth_DIF
InputData_Smth_DIF[][]=(InputData_Smth_DIF[p][q]<0)?(0):(InputData_Smth_DIF[p][q])

// Parsing response Averages to the position map
make /o/n=(nPositions,nROIs) SwooshMaps = NaN
make /o/n=(MapX_Max - MapX_Min+MapPxSize*2, MapY_Max-MapY_min,nRois) SwooshMaps_Px = 0
make /o/n=(MapX_Max - MapX_Min+MapPxSize*2, MapY_Max-MapY_min,nRois) SwooshMaps_Px0 = 0


variable nP_per_position = Spotdur_s/LineDuration

 
for (rr=0;rr<nROIs;rr+=1)

	for (pp=1;pp<nPositions-1;pp+=1)	 // make response per point array

		make /o/n=(nP_per_position) currentwave = InputData_Smth_DIF[p+pp*nP_per_position][rr]
		WaveStats/Q currentwave
	
		
		SwooshMaps[pp][rr] = V_SDev
	
		// find out average position of the swoop between points
		variable oldX = SwooshPositions[pp-1][1]
		variable oldY = SwooshPositions[pp-1][0]		
		variable newX = SwooshPositions[pp][1]
		variable newY = SwooshPositions[pp][0]
	
		variable midX = (newX-oldX)/2+oldX 
		variable midY = (newY-oldY)/2+oldY	
	
		// map that response to the pixel-map
		variable MapXPos1 = midX - MapX_Min 
		variable MapXPos2 = midX - MapX_Min  + 2*MapPxSize
		variable MapYPos1 = midY - MapY_Min  
		variable MapYPos2 = midY - MapY_Min  + 2*MapPxSize
	
		SwooshMaps_Px0[MapXPos1,MapXPos2][MapYPos1,MapYPos2][rr]+=1
		SwooshMaps_Px[MapXPos1,MapXPos2][MapYPos1,MapYPos2][rr]+=V_SDev
		
	endfor

	SwooshMaps_Px[][][rr]/=SwooshMaps_Px0[p][q][rr] // normalise maps by number of passes
	
endfor

// export handling
for (rr=0;rr<nROIs;rr+=1)
	string output_name2 = 'output_name1'+Num2Str(rr)
	make /o/n=(MapX_Max - MapX_Min+MapPxSize*2, MapY_Max-MapY_min) currentmap = SwooshMaps_Px[p][q][rr]
	setscale x,MapX_Min-MapPxSize,MapX_Max+MapPxSize,"mA" currentmap
	setscale y,MapY_Min-MapPxSize,MapY_Max+MapPxSize,"mA" currentmap
	duplicate /o Currentmap $output_name2
endfor
setscale x,MapX_Min-MapPxSize,MapX_Max+MapPxSize,"mA" SwooshMaps_Px
setscale y,MapY_Min-MapPxSize,MapY_Max+MapPxSize,"mA" SwooshMaps_Px

duplicate /o SwooshMaps_Px $output_name1
	
// display

if (Display_maps==1)
	for (rr=0;rr<nROIs;rr+=1)
		display /k=1
		ModifyGraph width=100
		ModifyGraph height={Aspect,0.6}
		string Displaymap = 'Output_name1'+Num2Str(rr)
		Appendimage $Displaymap
		•ModifyGraph mirror=0,fSize=8,axisEnab(left)={0.05,1},axisEnab(bottom)={0.05,1};DelayUpdate
		•Label left "\\Z10Y-Current (\U)";DelayUpdate
		•Label bottom "\\Z10ROI "+Num2Str(rr)+", X-Current (\U)"
		string mapname = "Swooshmaps"+Num2Str(Channel)+Num2Str(rr)
		ModifyImage $mapname ctab= {*,*,Rainbow,1}
		SetAxis/A/R left;DelayUpdate
		SetAxis/A/R bottom
	endfor
endif


// cleanup
//killwaves InputData, InputData_shift, InputData_smth, InputData_smth_DIF, SwooshMaps, Swooshmaps_Px, SwooshMaps_Px0, currentwave, Currentmap


end

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