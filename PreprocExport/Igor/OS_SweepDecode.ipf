#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_SweepDecode(display_stuff)
variable display_stuff  

wave SweepDecode
wave SkittlesWavelengths
wave Snippets0
wave SnippetsTimes0

variable nP_evaluatedelay = 0 // 25 == 50 ms
variable LineDuration_s = 0.002 // s
variable InstanceDuration_s = 0.4 // s

variable nLEDs = Dimsize(SweepDecode,0)
variable nWithinSweepRepeats = Dimsize(SweepDecode,1)
variable nInstances = nLEDs * nWithinSweepRepeats

variable nP_Snippet = Dimsize(Snippets0,0)
variable nLoops_Snippet = Dimsize(Snippets0,1)
variable nROIs = Dimsize(Snippets0,2)

variable nP_instance = InstanceDuration_s / LineDuration_s

variable nActualRepeats = nLoops_Snippet*nWithinSweepRepeats

variable rr, ROI, LED, ll

/////////// Break all into snippets
make /o/n=(nP_instance,nActualRepeats,nLEDs,nROIs) SweepSnippets = NaN
make /o/n=(nP_instance,nLEDs,nROIs) SweepAverages = 0

make /o/n=(nLEDs,nActualRepeats,nROIs) Sweep_Tuning_Snippets= 0
make /o/n=(nLEDs,nROIs) Sweep_Tuning_Ave = 0


for (ROI=0;ROI<nROIs;ROI+=1)
	for (LED=0;LED<nLEDs;LED+=1)
		for (ll=0;ll<nLoops_Snippet;ll+=1)
			for (rr=0;rr<nWithinSweepRepeats;rr+=1)

				variable currentTimePosition = SweepDecode[LED][rr]
				variable startpoint = currentTimePosition * nP_instance + nP_evaluatedelay
				
				SweepSnippets[][ll*nWithinSweepRepeats+rr][LED][ROI]=Snippets0[startpoint+p][ll][ROI]
				SweepAverages[][LED][ROI]+=Snippets0[startpoint+p][ll][ROI] / nActualRepeats
		
				make /o/n=(nP_instance) currentwave = SweepSnippets[p][ll*nWithinSweepRepeats+rr][LED][ROI]
				Wavestats/Q CurrentWave
				Sweep_Tuning_Snippets[LED][ll*nWithinSweepRepeats+rr][ROI]=V_SDev
	
			endfor
		endfor
		make /o/n=(nActualRepeats) currentwave = Sweep_Tuning_Snippets[LED][p][ROI]
		Wavestats/Q CurrentWave
		Sweep_Tuning_Ave[LED][ROI]=V_Avg

	endfor
endfor
Setscale /p x,0,LineDuration_s,"s" SweepSnippets, SweepAverages

// make one long sweep trace with snippets and averages
make /o/n=(nP_instance*nLEDs,nActualRepeats,nROIs) SweepSortedSnippets = NaN
make /o/n=(nP_instance*nLEDs,nROIs) SweepSortedAverages = 0

variable ii
for (ROI=0;ROI<nROIs;ROI+=1)
	for (LED=0;LED<nLEDs;LED+=1)
		for (ii=0;ii<nActualRepeats;ii+=1)
			startpoint = (nLEDs-1-LED) *  nP_instance
			SweepSortedSnippets[startpoint,startpoint+nP_instance][ii][ROI]=SweepSnippets[p-startpoint][ii][LED][ROI]
		endfor
		SweepSortedAverages[startpoint,startpoint+nP_instance][ROI]=SweepAverages[p-startpoint][LED][ROI]			
	endfor
endfor

Setscale /p x,0,LineDuration_s,"s" SweepSortedSnippets, SweepSortedAverages


///// DISPLAY

// make Stim wave
make /o/n=(nP_instance*nLEDs,nLEDs) SweepStimArray = 0
for (LED=0;LED<nLEDs;LED+=1)
	SweepStimArray[LED*nP_instance,(LED+0.5)*nP_instance][LED]=1
endfor
Setscale /p x,0,LineDuration_s,"s" SweepStimArray


if (display_stuff==1)

	// Trace display
	display /k=1
	make /o/n=(1) M_Colors
	Colortab2Wave Rainbow256
	string YAxisName
	string Tracename

	for (LED=0;LED<nLEDs;LED+=1)
		Appendtograph /l=StimY SweepStimArray[][LED]
	endfor
	ModifyGraph noLabel(StimY)=2,axThick(StimY)=0,axisEnab(StimY)={0.05,1};DelayUpdate
	ModifyGraph freePos(StimY)={0,kwFraction}


	for (ROI=0;ROI<nROIs;ROI+=1)
		YAxisName = "Y_ROI"+Num2Str(ROI)
		for (ii=0;ii<nActualRepeats;ii+=1)
			Appendtograph /l=$YAxisName SweepSortedSnippets[][ii][ROI]
		endfor
		ModifyGraph freePos($YAxisName)={0,kwFraction}
	endfor
	ModifyGraph rgb=(52224,52224,52224)
	ModifyGraph fSize=8,axisEnab(bottom)={0.05,1}
	
	for (LED=0;LED<nLEDs;LED+=1)
		Tracename = "SweepStimArray#"+Num2Str(LED)
		if (LED==0)
			Tracename = "SweepStimArray"
		endif
		variable colorposition = 256-(255 * (LED+1)/nLEDs)
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
		ModifyGraph mode($tracename)=7,hbFill($tracename)=5

	endfor

	for (ROI=0;ROI<nROIs;ROI+=1)
		YAxisName = "Y_ROI"+Num2Str(ROI)
		Tracename = "SweepSortedAverages#"+Num2Str(ROI)
		if (ROI==0)
			Tracename = "SweepSortedAverages"
		endif
		Appendtograph /l=$YAxisName SweepSortedAverages[][ROI]
		ModifyGraph lsize($TraceName)=1.5
		colorposition = 255 * (ROI+1)/nRois
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])

		variable plotfrom = (1-((ROI+1)/nRois))*0.8+0.2
		variable plotto = (1-(ROI/nRois))*0.8+0.2
		ModifyGraph axisEnab($YAxisName)={plotfrom,plotto}
		Label $YAxisName "\\Z10"+Num2Str(ROI)
		ModifyGraph noLabel($YAxisName)=1,axThick($YAxisName)=0
		ModifyGraph lblRot($YAxisName)=-90
	endfor
	
	// Evaluated display 
	
	display /k=1

	for (ROI=0;ROI<nROIs;ROI+=1)
		YAxisName = "Y_ROI"+Num2Str(ROI)
		for (ii=0;ii<nActualRepeats;ii+=1)
			Appendtograph /l=$YAxisName Sweep_Tuning_Snippets[][ii][ROI] vs SkittlesWavelengths
		endfor
		ModifyGraph freePos($YAxisName)={0,kwFraction}
	endfor
	ModifyGraph rgb=(52224,52224,52224)
	ModifyGraph fSize=8,axisEnab(bottom)={0.05,1}

	for (ROI=0;ROI<nROIs;ROI+=1)
		YAxisName = "Y_ROI"+Num2Str(ROI)
		Tracename = "Sweep_Tuning_Ave#"+Num2Str(ROI)
		if (ROI==0)
			Tracename = "Sweep_Tuning_Ave"
		endif
		Appendtograph /l=$YAxisName Sweep_Tuning_Ave[][ROI] vs SkittlesWavelengths
		ModifyGraph lsize($TraceName)=1.5
		colorposition = 255 * (ROI+1)/nRois
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])

		plotfrom = (1-((ROI+1)/nRois))*0.8+0.2
		plotto = (1-(ROI/nRois))*0.8+0.2
		ModifyGraph axisEnab($YAxisName)={plotfrom,plotto}
		Label $YAxisName "\\Z10"+Num2Str(ROI)
		ModifyGraph noLabel($YAxisName)=1,axThick($YAxisName)=0
		ModifyGraph lblRot($YAxisName)=-90
	endfor
	Label bottom "\\Z10Wavelength (nm)"


endif

print "Warning, the SnippetTime array is not yet implemented to trace times are not offset for scan Y position"

end

////////////////////////////////////////////////

function AddSweepStimTrace()

wave SweepDecode
wave Snippets0

variable LineDuration_s = 0.002 // s
variable InstanceDuration_s = 0.4 // s

variable nLEDs = Dimsize(SweepDecode,0)
variable nWithinSweepRepeats = Dimsize(SweepDecode,1)
variable nInstances = nLEDs * nWithinSweepRepeats

variable nP_Snippet = Dimsize(Snippets0,0)
variable nLoops_Snippet = Dimsize(Snippets0,1)
variable nROIs = Dimsize(Snippets0,2)

variable nP_instance = InstanceDuration_s / LineDuration_s

variable rr, ROI, LED, ll

/// Calculate original Stimulus Array for superposition to input avergage data
make /o/n=(nP_instance*nInstances, nLEDs) SweepStimArray_original = 0
for (LED=0;LED<nLEDs;LED+=1)
	for (rr=0;rr<nWithinSweepRepeats;rr+=1)
		variable Currentposition  = SweepDecode[LED][rr]
		SweepStimArray_original[Currentposition*nP_Instance,(Currentposition+0.5)*nP_Instance-1][LED]=1
	endfor
endfor
Setscale /p x,0,LineDuration_s,"s" SweepStimArray_original

//// DISPLAY


// Trace display
make /o/n=(1) M_Colors
Colortab2Wave Rainbow256
string Tracename

for (LED=0;LED<nLEDs;LED+=1)
	Appendtograph /l=StimY2 SweepStimArray_original[][LED]
	Tracename = "SweepStimArray_original#"+Num2Str(LED)
	if (LED==0)
		Tracename = "SweepStimArray_original"
	endif
	variable colorposition = (255 * (LED+1)/nLEDs)
	ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
	ModifyGraph mode($tracename)=7,hbFill($tracename)=5
	ReorderTraces AverageStimArtifact0,{$tracename}
endfor
ModifyGraph noLabel(StimY2)=2,axThick(StimY2)=0,axisEnab(StimY2)={0.05,1};DelayUpdate
ModifyGraph freePos(StimY2)={0,kwFraction}


end