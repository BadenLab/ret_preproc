#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_CallManualRoi()

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

// flags from "OS_Parameters"
variable X_cut = OS_Parameters[%LightArtifact_cut]
variable LineDuration = OS_Parameters[%LineDuration]

// data handling
wave wParamsNum // Reads data-header
string input_name = "wDataCh"+Num2Str(Channel)+"_detrended"
duplicate /o $input_name InputData
variable nX = DimSize(InputData,0)
variable nY = DimSize(InputData,1)
variable nF = DimSize(InputData,2)
variable Framerate = 1/(nY * LineDuration) // Hz 
variable Total_time = (nF * nX ) * LineDuration
print "Recorded ", total_time, "s @", framerate, "Hz"
variable xx,yy,ff // initialise counters

// calculate Pixel / ROI sizes in microns
variable zoom = wParamsNum(30) // extract zoom
variable px_Size = (0.65/zoom * 110)/nX // microns
print "Pixel Size:", round(px_size*100)/100," microns"

// make SD average
make /o/n=(nX,nY) Stack_SD = 0 // Avg projection of InputData
make /o/n=(nX,nY) ROIs = 1 // empty ROI wave
make /o/n=(nF) currentwave = 0
for (xx=X_cut;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		Multithread currentwave[]=InputData[xx][yy][p] // get trace from "reference pixel"
		Wavestats/Q currentwave
		Stack_SD[xx][yy]=V_SDev
	endfor
endfor
setscale /p x,-nX/2*px_Size,px_Size,"µm" ROIs, Stack_SD
setscale /p y,-nY/2*px_Size,px_Size,"µm"  ROIs, Stack_SD

// display SD wave
Display /k=1 
Appendimage Stack_SD
Appendimage ROIs
ModifyImage ROIs explicit=1,eval={-1,65535,0,0}
ModifyGraph height={Aspect,nY/nX}
WMCreateImageROIPanel() // calls SARFIA Roi generator - if follow that through it gives a wave 

// cleanup
killwaves currentwave,InputData

end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_ApplyManualRoi()

// data handling
wave M_ROIMask
variable nX = Dimsize(M_ROIMask,0)
variable nY = Dimsize(M_ROIMask,1)
make /o/n=(nX,nY) ROIs = 1 // empty ROI wave

variable xx,yy,rr

// calculate Pixel size in microns to scale ROIs
wave wParamsNum // Reads data-header
variable zoom = wParamsNum(30) // extract zoom
variable px_Size = (0.65/zoom * 110)/nX // microns
setscale /p x,-nX/2*px_Size,px_Size,"µm" ROIs
setscale /p y,-nY/2*px_Size,px_Size,"µm"  ROIs

// create proper ROI Mask from M_ROIMask
duplicate /o M_ROIMask ROIbw_sub // make a lookup wave
duplicate /o ROIs ROIs_add
variable nRois = 0
for (xx=0;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)
		if (ROIbw_sub[xx][yy]==0) // ROI are coded as 0, space is coded as 1...
			// here it found a seed
			nRois+=1
			variable RoiValue = (nRois)*(-1)
			ROIs_add[xx][yy]=RoiValue // add it to the ROI wave
			ROIs[xx][yy]=RoiValue // add it to the ROI wave
			do // flood fill that ROI
				if (nY>1)
					Imagestats/Q ROIs
				elseif (nY==1)
					Wavestats/Q ROIs // otherwise gets error from imagestats
				endif
				
				variable ROI_average_before = V_Avg
				// flood fill: note this will always work but if a ROI touches the edge it will report an out of range error.just live with it, otherwise 
				// would need endless set of if clauses
				
				if (nY>1)
					Multithread ROIs_add[0,nX-1][0,nY-1]=((ROIbw_sub[p][q]==0) && ((ROIs[p+1][q]==RoiValue)||(ROIs[p-1][q]==RoiValue)||(ROIs[p][q+1]==RoiValue)||(ROIs[p][q-1]==RoiValue)))?(RoiValue):(ROIs[p][q]) // flood fill
				elseif (nY==1) // no Y dimension in linescan)
					Multithread ROIs_add[0,nX-1][0]=((ROIbw_sub[p][q]==0) && ((ROIs[p+1][q]==RoiValue)||(ROIs[p-1][q]==RoiValue)))?(RoiValue):(ROIs[p][q]) // flood fill
				endif


				DoUpdate
				ROIs = ROIs_add
				if (nY>1)
					Imagestats/Q ROIs
				elseif (nY==1)
					Wavestats/Q ROIs // otherwise gets error from imagestats
				endif
				
				variable ROI_average_after = V_Avg
				if (ROI_average_before==ROI_average_after) // leaves fill if no more change
					break
				endif			
			while(1)
			ROIbw_sub[][]=(ROIs[p][q]==RoiValue)?(1):(ROIbw_sub[p][q]) // kill that ROI from lookup wave
		endif
	endfor
endfor

// colour in the ROIs
make /o/n=(1) M_Colors
Colortab2Wave Rainbow256
for (rr=0;rr<nRois;rr+=1)
	variable colorposition = 255 * (rr+1)/nRois
	ModifyImage ROIs explicit=1,eval={-rr-1,M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2]}
endfor



// cleanup
killwaves ROIbw_sub, ROIs_add,M_Colors

print RoiValue*(-1) ,"ROIs generated manually"


end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_monoPixelApply()

// data handling
wave M_ROIMask
variable nX = Dimsize(M_ROIMask,0)
variable nY = Dimsize(M_ROIMask,1)
make /o/n=(nX,nY) ROIs = 1 // empty ROI wave

variable xx,yy,rr

// calculate Pixel size in microns to scale ROIs
wave wParamsNum // Reads data-header
variable zoom = wParamsNum(30) // extract zoom
variable px_Size = (0.65/zoom * 110)/nX // microns
setscale /p x,-nX/2*px_Size,px_Size,"µm" ROIs
setscale /p y,-nY/2*px_Size,px_Size,"µm"  ROIs

// make each pixel == 1 ROI
variable nRois = 0

for (xx=0;xx<nX;xx+=1)
	for (yy=0;yy<nY;yy+=1)

		if (M_ROIMask[xx][yy]==0) // if its part of the marked region
			ROIs[xx][yy]=-nRois-1
			nRois+=1
		endif
	endfor
endfor

print nRois, "pixels selected as individual ROIs"
// colour in the ROIs
make /o/n=(1) M_Colors
Colortab2Wave Rainbow256
for (rr=0;rr<nRois;rr+=1)
	variable colorposition = 255 * (rr+1)/nRois
	ModifyImage ROIs explicit=1,eval={-rr-1,M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2]}
endfor



// cleanup
killwaves M_Colors



end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_CloneSarfiaRoi()

// 1 // check for existing Sarfia Mask from Laplace Operator thing
if (waveexists($"ROIs")==0) // KF 20150310; changed name of preexisting roi mask from MTROIWave to ROIs
	print "Warning: ROIs does not exist - doing nothing..."
else	
	
	wave ROIs
		
	wave OS_Parameters
	// flags from "OS_Parameters"
	variable X_cut = OS_Parameters[%LightArtifact_cut]
	variable LineDuration = OS_Parameters[%LineDuration]
	variable Channel = OS_Parameters[%Data_Channel]
	
	// data handling
	wave wParamsNum // Reads data-header
	string input_name = "wDataCh"+Num2Str(Channel)+"_detrended"
	duplicate /o $input_name InputData
	variable nX = DimSize(InputData,0)
	variable nY = DimSize(InputData,1)
	variable nF = DimSize(InputData,2)
	variable Framerate = 1/(nY * LineDuration) // Hz 
	variable Total_time = (nF * nX ) * LineDuration
	print "Recorded ", total_time, "s @", framerate, "Hz"
	
	variable xx,yy,rr
	
	// calculate Pixel / ROI sizes in microns
	variable zoom = wParamsNum(30) // extract zoom
	variable px_Size = (0.65/zoom * 110)/nX // microns
	print "Pixel Size:", round(px_size*100)/100," microns"
	
	// make SD average
	make /o/n=(nX,nY) Stack_SD = 0 // Avg projection of InputData
	make /o/n=(nX,nY) ROIs_corrected = 1 // empty ROI wave
	make /o/n=(nF) currentwave = 0
	for (xx=X_cut;xx<nX;xx+=1)
		for (yy=0;yy<nY;yy+=1)
			Multithread currentwave[]=InputData[xx][yy][p] // get trace from "reference pixel"
			Wavestats/Q currentwave
			Stack_SD[xx][yy]=V_SDev
		endfor
	endfor
	
	// take MTROIWave and work out if it's dimensions are cropped in X axis (e.g. to get rid of light artifact before using the Laplace)
	Variable nX_ROIs = Dimsize(ROIs,0)
	
	// Project MTROIwave onto ROIs with correct x offset
	ROIs_corrected[nX-nX_ROIs,nX-1][]=ROIs[p-(nX-nX_ROIs)][q]
	duplicate/o ROIs_corrected, ROIs
	killwaves ROIs_corrected
	
	setscale /p x,-nX/2*px_Size,px_Size,"µm" Stack_SD,ROIs
	setscale /p y,-nY/2*px_Size,px_Size,"µm" Stack_SD, ROIs
	
	// display SD wave
	Display /k=1
	Appendimage Stack_SD
	Appendimage ROIs
	ModifyImage ROIs explicit=1,eval={-1,65535,0,0}
	ModifyGraph height={Aspect,nY/nX}
	
	// colour in the ROIs
	make /o/n=(1) M_Colors
	Colortab2Wave Rainbow256
	variable nRois = Wavemin(ROIs)*(-1)
	for (rr=0;rr<nRois;rr+=1)
		variable colorposition = 255 * (rr+1)/nRois
		ModifyImage ROIs explicit=1,eval={-rr-1,M_Colors[colorposition][0],M_Colors[colorposition][1],M_Colors[colorposition][2]}
	endfor
	
	print "ROIs generated from SARFIA Mask"
		
	// cleanup
	killwaves currentwave,InputData

endif

end





end