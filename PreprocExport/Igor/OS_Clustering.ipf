#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_Clustering()

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
variable Display_averages = OS_Parameters[%Display_Stuff]
variable LineDuration = OS_Parameters[%LineDuration]
variable nClasses_seed = OS_Parameters[%Clustering_nClasses]
variable nSD_plot = OS_Parameters[%Clustering_SDplot]

// data handling
string input_name = "Averages"+Num2Str(Channel)
duplicate /o $input_name InputData

variable nP = DimSize(InputData,0)
variable nRois = DimSize(InputData,1)

wave AverageStimArtifact0 // generated by "Traces and Triggers"
wave ROIs
wave Stack_SD
wave qualitycriterion


string output_name1 = "Cluster_Means"+Num2Str(Channel)
string output_name2 = "Cluster_Allocations"+Num2Str(Channel)

variable pp,rr,cc

//// Kick out bad cells
//Variable QualityThreshold = 0.2
//
//Variable BadROIs = 0
//for (rr=0;rr<nRois;rr+=1)
//	if (QualityCriterion[rr]<QualityThreshold)
//		InputData[][rr]=0
//		BadROIs+=1
//	endif
//endfor
//Print "Clustered ",nROIs-BadROIs, "/", nROIs, "ROIs"

// Clustering
make /o/n=1 M_KMClasses
make /o/n=(nROIs) W_KMMembers
Kmeans /NCLS=(nClasses_seed) /OUT=2 InputData
// now M_KMClasses has the cluster means (classes) and W_KMMembers has each ROI's class assignment
variable nClasses = Dimsize(M_KMClasses,1)
setscale /p x,0,LineDuration,"s" M_KMClasses
M_KMClasses[nP-1][] = NaN


// export handling
duplicate /o M_KMClasses $output_name1
duplicate /o W_KMMembers $output_name2


	
// display

if (Display_averages==1)
	display /k=1
	
	make /o/n=(1) M_Colors
	Colortab2Wave Rainbow256

	Appendtograph /l=StimY AverageStimArtifact0
	ModifyGraph fSize=8,noLabel(StimY)=2,axThick(StimY)=0,lblPos(StimY)=47;DelayUpdate
	ModifyGraph axisEnab(StimY)={0.05,0.15},freePos(StimY)={0,kwFraction}
	ModifyGraph rgb(AverageStimArtifact0)=(0,0,0)
	
	string tracename
	variable rr_counter =0
	for (cc=0;cc<nClasses;cc+=1)
		string YAxisName = "YClass"+Num2Str(cc)
		// add the individual members of a class
		for (rr=0;rr<nRois;rr+=1)
			if (W_KMMembers[rr]==cc)
				string tracename2 = input_name+"#"+Num2Str(rr_counter)
				Appendtograph /l=$YAxisName $input_name[][rr] 
				if (rr_counter==0)
					tracename2 = input_name
				endif
				ModifyGraph rgb($tracename2)=(52224,52224,52224)
				rr_counter+=1
			endif
		endfor
		
		// Add the class
		tracename = output_name1+"#"+Num2Str(cc)
		if (cc==0)
			tracename = output_name1
		endif
		Appendtograph /l=$YAxisName $output_name1[][cc]
		variable colorposition = 255 * (cc+1)/nClasses
		ModifyGraph rgb($tracename)=(M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2])
		ModifyGraph lsize($tracename)=1.5

		variable plotfrom = (1-((cc+1)/nClasses))*0.8+0.2
		variable plotto = (1-(cc/nClasses))*0.8+0.2
		ModifyGraph fSize($YAxisName)=8,axisEnab($YAxisName)={plotfrom,plotto};DelayUpdate
		ModifyGraph freePos($YAxisName)={0,kwFraction};DelayUpdate
		Label $YAxisName "\\Z10"+Num2Str(cc)
		ModifyGraph noLabel($YAxisName)=1,axThick($YAxisName)=0;DelayUpdate
		ModifyGraph lblRot($YAxisName)=-90
		
		SetAxis $YAxisName -nSD_plot,nSD_plot
	endfor
	ModifyGraph fSize=8
	
	�ModifyGraph axisEnab(bottom)={0.05,1};DelayUpdate
	�Label bottom "\\Z10Time (\\U)"

	// RoiMap plotting
	
	display /k=1 
	Appendimage Stack_SD
	Appendimage ROIs
	for (rr=0;rr<nRois;rr+=1)
		colorposition = 255 * (W_KMMembers[rr]+1)/nClasses
		ModifyImage ROIs explicit=1,eval={-rr,M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2]}
	endfor
	
	

endif


// cleanup
killwaves InputData, M_KMClasses, W_KMMembers


end