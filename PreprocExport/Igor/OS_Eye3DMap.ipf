#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_Eye3DMap(display_stuff)
variable display_stuff

// MANUAL KEY PARAMETERS

variable AnimalAngle_deg = 0
variable EyeLeftRight = -1 // left eye 1, right eye -1
wave wParamsNum
// change the "+0" bit if the sutter was moved relative to eye origin
variable Origin_SutterAbsX =  wParamsNum[FindDimLabel(wParamsNum,0,"XCoord_um")] +0
variable Origin_SutterAbsY =  wParamsNum[FindDimLabel(wParamsNum,0,"YCoord_um")] + 0
variable Origin_SutterAbsZ =  wParamsNum[FindDimLabel(wParamsNum,0,"ZCoord_um")] + 0

variable display_3D = 1

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ROTATION AND POSITION NOTES, BASED ON SETUP 1
// Up in image (Y) increses Sutter X position
// Right in image (X) increases Sutter Y position
// Down in focus (Z) increases Sutter Z position
// positive angle rotation goes anticlockwise in image
// Xdrag offset in image is positive Voffset in X coordinate
// Ydrag offset in image is positive Voffset in Y coordinate
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// get/calculate key parameters from header info (automatic)

variable VtoDistanceFactor = 82 // 82 microns per Volt under x20 Objective
variable IPL_to_micron_scale = 30 // how thick is the IPL
variable AnimalAngle = AnimalAngle_deg / 360  * pi 

wave positions
wave ROIs
wave CoM

variable Scan_SutterAbsX =  wParamsNum[FindDimLabel(wParamsNum,0,"XCoord_um")]
variable Scan_SutterAbsY =  wParamsNum[FindDimLabel(wParamsNum,0,"YCoord_um")]
variable Scan_SutterAbsZ =  wParamsNum[FindDimLabel(wParamsNum,0,"ZCoord_um")]

variable SutterOffsetX = round(Scan_SutterAbsX - Origin_SutterAbsX)
variable SutterOffsetY = round(Scan_SutterAbsY - Origin_SutterAbsY)
variable SutterOffsetZ = round(Scan_SutterAbsZ - Origin_SutterAbsZ)

variable Zoom =  wParamsNum[FindDimLabel(wParamsNum,0,"Zoom")]
variable ScanAngle_deg = wParamsNum[FindDimLabel(wParamsNum,0,"Angle_deg")]
variable ScanAngle = (ScanAngle_deg / 360)  * (2*pi)

variable X_UserOffsetV =  wParamsNum[FindDimLabel(wParamsNum,0,"User_XOffset_V")]
variable Y_UserOffsetV =  wParamsNum[FindDimLabel(wParamsNum,0,"User_YOffset_V")]

variable xPixelsInd = FindDimLabel(wParamsNum,0,"User_dxPix" )
variable yPixelsInd = FindDimLabel(wParamsNum,0,"User_dyPix" )
variable realPixDurInd = FindDimLabel(wParamsNum,0,"RealPixDur" )
variable lineDur = (wParamsNum[xPixelsInd] *  wParamsNum[realPixDurInd]) * 10^-6

variable nROIs = Dimsize(positions,0)
variable nX = Dimsize(ROIs,0)
variable nLines = Dimsize(ROIs,1)
variable px_to_microns = 0.65 / Zoom * 110  / nX

//////////////////////////////////////

// make standard Eye semiSphere - for 3D display (optional)
variable sphereScale = 150
variable nTemplatePoints = 1000
make /o/n=(nTemplatePoints,3) Eye3DTemplate_xyz = NaN
variable offset = 1 / nTemplatePoints
variable increment = pi * (3 - sqrt(5))
variable phi=0
variable pp
for (pp=0;pp<nTemplatePoints;pp+=1)
	variable sphereZ = (((pp * offset) -1) + (offset / 2))*-1
	variable radius = (1-sphereZ^2)^0.5
	phi+=increment
	variable sphereX = cos(phi) * radius
	variable sphereY = sin(phi) * radius
	Eye3DTemplate_xyz[pp][0]=sphereX * sphereScale
	Eye3DTemplate_xyz[pp][1]=sphereY * sphereScale
	Eye3DTemplate_xyz[pp][2]=sphereZ * sphereScale		
endfor
make /o/n=1 M_interpolatedImage = NaN
imageinterpolate /s={-sphereScale,1,sphereScale,-sphereScale,1,sphereScale} Voronoi Eye3DTemplate_xyz
duplicate /o M_interpolatedImage Eye3DTemplate_surface
killwaves M_interpolatedImage

/// for each ROI get both in cartesian coordinates (xyz: Positions_3D_xyz) and vector (angle + IPL depth: Positions_2D)

make /o/n=(nRois,3) Positions_3D_xyz = NaN // X, Y, Z
make /o/n=(nRois,3) Positions_2D = NaN // angle and IPL depth

variable Scan_X = Scan_SutterAbsY - Origin_SutterAbsY + X_UserOffsetV*VtoDistanceFactor // note that X and Y is swapped as the Sutter is flipped
variable Scan_Y = Scan_SutterAbsX - Origin_SutterAbsX + Y_UserOffsetV*VtoDistanceFactor
variable Scan_Z = Scan_SutterAbsZ - Origin_SutterAbsZ  // Z lens not implemented yet

variable rr
for (rr=0;rr<nRois;rr+=1)

	variable ROIinScan_X = ((nX/2 - CoM[rr][0] ) * px_to_microns ) *-1 // X inverted
	variable ROIinScan_Y = ((nLines/2 - CoM[rr][1] ) * px_to_microns) * -1 // Y inverted as ScanM plots them upside down
	
	variable ROIinScan_X_rotated = (ROIinScan_X) * cos(ScanAngle) - (ROIinScan_Y) * sin(ScanAngle)
	variable ROIinScan_Y_rotated = (ROIinScan_X) * sin(ScanAngle) + (ROIinScan_Y) * cos(ScanAngle)
	
	variable ROIinScan_X_rotated_offset = ROIinScan_X_rotated + Scan_X
	variable ROIinScan_Y_rotated_offset = ROIinScan_Y_rotated + Scan_Y

	variable ROI_X = ((ROIinScan_X_rotated_offset) * cos(-AnimalAngle) - (ROIinScan_Y_rotated_offset) * sin(-AnimalAngle) ) * EyeLeftRight
	variable ROI_Y = (ROIinScan_X_rotated_offset) * sin(-AnimalAngle) + (ROIinScan_Y_rotated_offset) * cos(-AnimalAngle)
	variable ROI_Z = Scan_Z // z lens not implemented yet		


	variable xyangle = atan2(ROI_Y,ROI_X)
	variable xyradius = (Scan_X^2+Scan_Y^2)^0.5 // note that only the xy angle is based on the ROIS, the rest is scan. this assumed that the scan is always horizontal (no z lens)
	variable zangle = atan2(xyradius,Scan_Z)

	Positions_2D[rr][0] = xyangle // angle
	Positions_2D[rr][1] = positions[rr] // IPL depth

	Positions_3D_xyz[rr][0]=cos(xyangle) * sin(zangle) * (xyradius + (positions[rr]/100-0.5) * IPL_to_micron_scale)
	Positions_3D_xyz[rr][1]=sin(xyangle) * sin(zangle) * (xyradius + (positions[rr]/100-0.5) * IPL_to_micron_scale)
	Positions_3D_xyz[rr][2]=Scan_Z

endfor


// DISPLAY

if (display_stuff==1)	

	if (display_3D==1)
		string GizmoCommand = "NewGizmo /k=1"
		Execute/Q GizmoCommand
		GizmoCommand = "AppendToGizmo nextSurface=Eye3DTemplate_surface"
		Execute/Q GizmoCommand	
		GizmoCommand = "AppendToGizmo DefaultScatter=Positions_3D_xyz"
		Execute/Q GizmoCommand
		GizmoCommand = "ModifyGizmo setOuterBox={-150,150,-150,150,-150,150}, euler = {180, 0, 0 }, aspectRatio = 1"
		Execute/Q GizmoCommand
	
		GizmoCommand = "ModifyGizmo compile"
		Execute/Q GizmoCommand
	endif
	//
	Display /k=1 Positions_3D_xyz[][1] vs Positions_3D_xyz[][0]
	•ModifyGraph zero=2,fSize=8,noLabel=1,axisEnab(left)={0.05,1};DelayUpdate
	•ModifyGraph axisEnab(bottom)={0.05,1};DelayUpdate
	•SetAxis left -150,150;DelayUpdate
	•SetAxis bottom -150,150
	•ModifyGraph mode=2,rgb=(0,0,0)
	•ModifyGraph width={Aspect,1}
	ModifyGraph nticks=3,noLabel=0;DelayUpdate
	Label left "\\Z10µm";DelayUpdate
	Label bottom "\\Z10µm"

	ShowTools/A 
	SetDrawEnv xcoord= bottom, ycoord= left,fsize= 10
	DrawText 120,-130+200,"V"
	SetDrawEnv xcoord= bottom, ycoord= left,fsize= 10
	DrawText 120,-70+200,"D"
	SetDrawEnv xcoord= bottom, ycoord= left,fsize= 10
	DrawText 90,-100+200,"N"
	SetDrawEnv xcoord= bottom, ycoord= left,fsize= 10
	DrawText 150,-100+200,"T"	
	
	SetDrawEnv xcoord= rel, ycoord= rel,fsize= 10
	DrawText 0.02,0.98,"Z-offset: "+Num2Str(Scan_Z)+"µm"
	HideTools/A 

endif


end


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This function goes through all scan folders that are inside the folder "root:scans:" and takes the 3D positions as well as kernelList_SD thing
// to make big arrays from all scans in one

function OS_Eye3D_Collect(display_stuff)
variable display_stuff

SetDataFolder root:Scans

variable nScanFolders = countObjects(":",4)
make /o/n=(nScanFolders) ROIs_per_scan = NaN
make /o/n=(nScanFolders+1) ROIs_per_scan_cumul = 0

String ItemAddress
String ItemTargetName

variable ff
for (ff=0;ff<nScanFolders;ff+=1)
	String DF_Scan = GetIndexedObjNameDFR(GetDataFolderDFR(),4,ff)
	String DF_Scan_Full = "root:Scans:'"+'DF_Scan'+"':"
	
	// grab all interesting arrays
	ItemAddress = DF_Scan_Full+"Kernels_all"	// RawKernels
	ItemTargetName = "Kernels_all"+Num2Str(ff)	
	duplicate /o $ItemAddress $ItemTargetName

	ItemAddress = DF_Scan_Full+"Positions_3D_xyz"	
	ItemTargetName = "Positions_3D_xyz"+Num2Str(ff)	
	duplicate /o $ItemAddress $ItemTargetName
	
	ItemAddress = DF_Scan_Full+"Positions_2D"	
	ItemTargetName = "Positions_2D"+Num2Str(ff)	
	duplicate /o $ItemAddress $ItemTargetName
	
	ItemAddress = DF_Scan_Full+"KernelList_SD"	// For example here we are taking the array "kernelList_SD" and copying it over into root:scans"
	ItemTargetName = "KernelList_SD"+Num2Str(ff)	
	duplicate /o $ItemAddress $ItemTargetName
	
	ItemAddress = DF_Scan_Full+"KernelList_FFT"	// FFT
	ItemTargetName = "KernelList_FFT"+Num2Str(ff)	
	duplicate /o $ItemAddress $ItemTargetName
	
	ItemAddress = DF_Scan_Full+"KernelList_BiPhasic"	// BiPhasic
	ItemTargetName = "KernelList_BiPhasic"+Num2Str(ff)	
	duplicate /o $ItemAddress $ItemTargetName
	
	// count ROIs
	ROIs_per_scan[ff]=Dimsize($ItemTargetName,0)
	ROIs_per_scan_cumul[ff+1,nScanFolders]+=ROIs_per_scan[ff]
endfor
Wavestats/Q ROIs_per_scan
variable nROIs = V_Sum
variable nPoints_Kernel = 649 // hard coded kernel duration
print "nScans: ", nScanfolders
print "nROIs: ", nROIs

// Glue arrays together

make /o/n=(nROIs,3) Positions_3D_xyz = NaN
make /o/n=(nROIs,3) Positions_2D = NaN
make /o/n=(nROIs,4) KernelList_SD = NaN // declare the array that you picked above
make /o/n=(nROIs,4) KernelList_FFT = NaN // declare the array that you picked above
make /o/n=(nROIs,4) KernelList_BiPhasic = NaN // declare the array that you picked above
make /o/n=(nPoints_Kernel,4,nROIs) Kernels_all = NaN // declare the array that you picked above

for (ff=0;ff<nScanFolders;ff+=1)
	ItemTargetName = "Positions_3D_xyz"+Num2Str(ff)	
	duplicate /o $ItemTargetName Currentwave
	killwaves $ItemTargetName
	Positions_3D_xyz[ROIs_per_scan_cumul[ff],ROIs_per_scan_cumul[ff+1]-1][]=Currentwave[p-ROIs_per_scan_cumul[ff]][q]
	
	ItemTargetName = "Positions_2D"+Num2Str(ff)	
	duplicate /o $ItemTargetName Currentwave
	killwaves $ItemTargetName
	Positions_2D[ROIs_per_scan_cumul[ff],ROIs_per_scan_cumul[ff+1]-1][]=Currentwave[p-ROIs_per_scan_cumul[ff]][q]

	ItemTargetName = "KernelList_SD"+Num2Str(ff)	// fill array
	duplicate /o $ItemTargetName Currentwave
	killwaves $ItemTargetName
	KernelList_SD[ROIs_per_scan_cumul[ff],ROIs_per_scan_cumul[ff+1]-1][]=Currentwave[p-ROIs_per_scan_cumul[ff]][q]
	
	ItemTargetName = "KernelList_FFT"+Num2Str(ff)	// fill array
	duplicate /o $ItemTargetName Currentwave
	killwaves $ItemTargetName
	KernelList_FFT[ROIs_per_scan_cumul[ff],ROIs_per_scan_cumul[ff+1]-1][]=Currentwave[p-ROIs_per_scan_cumul[ff]][q]

	ItemTargetName = "KernelList_BiPhasic"+Num2Str(ff)	// fill array
	duplicate /o $ItemTargetName Currentwave
	killwaves $ItemTargetName
	KernelList_BiPhasic[ROIs_per_scan_cumul[ff],ROIs_per_scan_cumul[ff+1]-1][]=Currentwave[p-ROIs_per_scan_cumul[ff]][q]
	
	ItemTargetName = "Kernels_all"+Num2Str(ff)	// fill array
	duplicate /o $ItemTargetName Currentwave
	killwaves $ItemTargetName
	Kernels_all[][][ROIs_per_scan_cumul[ff],ROIs_per_scan_cumul[ff+1]-1]=Currentwave[p][q][r-ROIs_per_scan_cumul[ff]]
endfor

// display
if (display_stuff==1)
	display /k=1
	Appendtograph /l=R_Y positions_2D[][1] vs positions_2D[][0] // plot each ROI as a dot at its correct coordinate
	Appendtograph /l=G_Y positions_2D[][1] vs positions_2D[][0]
	Appendtograph /l=B_Y positions_2D[][1] vs positions_2D[][0]
	Appendtograph /l=U_Y positions_2D[][1] vs positions_2D[][0]	

	ModifyGraph zColor(Positions_2D)={KernelList_SD[*][0],-10,10,Red,0} // colour it in according to the picked array
	ModifyGraph zColor(Positions_2D#1)={KernelList_SD[*][1],-10,10,Green,0}
	ModifyGraph zColor(Positions_2D#2)={KernelList_SD[*][2],-10,10,Blue,0}
	ModifyGraph zColor(Positions_2D#3)={KernelList_SD[*][3],-10,10,Magenta,0}

	ModifyGraph mode=3,marker=16,msize=1			
	
	ModifyGraph nticks(R_Y)=2,nticks(G_Y)=2,nticks(B_Y)=2,nticks(U_Y)=2,fSize=8;DelayUpdate
	•ModifyGraph lblPos(R_Y)=47,lblPos(G_Y)=47,lblPos(B_Y)=47,lblPos(U_Y)=47;DelayUpdate
	•ModifyGraph axisEnab(R_Y)={0.8,1},axisEnab(bottom)={0.05,1};DelayUpdate
	•ModifyGraph axisEnab(G_Y)={0.55,0.75},axisEnab(B_Y)={0.3,0.5};DelayUpdate
	•ModifyGraph axisEnab(U_Y)={0.05,0.25},freePos(R_Y)={0,kwFraction};DelayUpdate
	•ModifyGraph freePos(G_Y)={0,kwFraction},freePos(B_Y)={0,kwFraction};DelayUpdate
	•ModifyGraph freePos(U_Y)={0,kwFraction};DelayUpdate
	•Label R_Y "\\Z10IPL (%)";DelayUpdate
	•Label bottom "\\Z10Eye position (rad.)";DelayUpdate
	•Label G_Y "\\Z10IPL (%)";DelayUpdate
	•Label B_Y "\\Z10IPL (%)";DelayUpdate
	•Label U_Y "\\Z10IPL (%)";DelayUpdate
	•SetAxis bottom -3.141,3.141
	
	display /k=1 Positions_3D_xyz[][1] vs Positions_3D_xyz[][0]
	•ModifyGraph zero=2,fSize=8,axisEnab(left)={0.05,1},axisEnab(bottom)={0.05,1};DelayUpdate
	•SetAxis left -150,150;DelayUpdate
	•SetAxis bottom -150,150
	•ModifyGraph mode=2,rgb=(0,0,0)
	•ModifyGraph width={Aspect,1}
	

endif

// cleanup
killwaves Currentwave

end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function OS_Eye3D_Analyse(display_stuff,iterateparameter1,iterateparameter2,iterateamplrange)
variable display_stuff,iterateparameter1,iterateparameter2,iterateamplrange
		
variable SD_threshold =4// consider only kernels above this SD
variable nBins_angle = 12 // 360 deg divided into how many sectors?
variable nBins_depth = 100 // 100% IPL depth divdied into how many positions
variable AngleSmooth = 60 // 60 // in degrees
variable DepthSmooth =6//6//5 // in % IPL


wave CI // get clusterIndex to pick the same cells as used in clustering
wave Positions_2D
wave KernelList_SD


//variable IPLE1_Ampl = 3// % IPL expansion at middle of range ("down")
//variable IPLE1_position = 0.86*(2*pi)
//variable IPLE1_bottom = 0// % IPL at bottom that ignores the expansion

//iterateparameter1 = 0.65
//iterateparameter2 = 5
//iterateamplrange = 20

variable IPLE_edgecut = 10

variable IPLE1_Ampl =5//iterateparameter2// % IPL expansion of bottom layers
variable IPLE1_position = 0.95 * (2*pi)// iterateparameter1*2*pi//iterateparameter1*(2*pi)

variable IPLE2_Ampl =5//iterateparameter2// % IPL expansion of top layers
variable IPLE2_position = 1 * (2*pi) // iterateparameter1*(2*pi)

variable IPLE_Ycut = 45// % IPL at bottom that ignores the compression

//IPLE1_Ampl = 0
//IPLE2_Ampl = 0

////////////////////////////////////////

//calculate expansion table 1 (top warp)
make /o/n=(nBins_angle,100) ExpansionTable1 = 0
Setscale x,0+IPLE1_position,2*pi+IPLE1_position,"" ExpansionTable1
ExpansionTable1[][0,IPLE_Ycut-1] = (y + IPLE1_Ampl * (Sin(x)+1)/2 ) * ((IPLE_Ycut)/(IPLE_Ycut+IPLE1_Ampl))

make /o/n=(nBins_angle,100) ExpansionTable2 = 0
Setscale x,0+IPLE2_position,2*pi+IPLE2_position,"" ExpansionTable2
ExpansionTable2[][IPLE_Ycut,99] = ((y-IPLE_Ycut) + (IPLE2_Ampl * (Sin(x)+1)/2 )) * ((100-IPLE_Ycut)/(IPLE2_Ampl+100-IPLE_Ycut)) +  IPLE_Ycut

make /o/n=(nBins_angle,100) ExpansionTable = 0
Setscale x,0,2*pi,"" ExpansionTable
ExpansionTable[][]=ExpansionTable1[p][q]+ExpansionTable2[p][q]


////////////////////////////////////////

variable offsetScale = (100/nBins_depth) / 30

AngleSmooth*=nBins_angle / 360
DepthSmooth*=nBins_depth / 100


variable nROIs = Dimsize(KernelList_SD,0)
variable nLEDs = 4

make /o/n=(nBins_angle,nBins_depth,nLEDs*2) Kernel_SDHists = 0
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_summed = 0
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_Rsummed = 0
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_Gsummed = 0
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_Bsummed = 0
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_Usummed = 0
make /o/n=(nBins_depth) Kernel_SDHists_summed2 = 0
make /o/n=(nBins_depth) Kernel_SDHists_Rsummed2 = 0
make /o/n=(nBins_depth) Kernel_SDHists_Gsummed2 = 0
make /o/n=(nBins_depth) Kernel_SDHists_Bsummed2 = 0
make /o/n=(nBins_depth) Kernel_SDHists_Usummed2 = 0

make /o/n=(nBins_angle,nBins_depth,nLEDs*2) Kernel_SDCounter = 1
make /o/n=(nBins_angle) Kernel_SDCounter2 = 1
make /o/n=(nBins_angle) Kernel_SDCounter3 = 1
make /o/n=(nBins_depth) Kernel_SDCounter_all = 1

variable rr,ll,dd
for (rr=0;rr<nROIs;rr+=1)
	variable CurrentAngle_bin = ((Positions_2D[rr][0] / (2*pi)) + 0.5) * nBins_angle

	variable CurrentIPLDepth = Positions_2D[rr][1]
	variable CurrentIPLDepth_corrected=ExpansionTable[CurrentAngle_bin][CurrentIPLDepth]

	CurrentIPLDepth_corrected/=100
	CurrentIPLDepth_corrected*=nBins_depth		
	
	for (ll=0;ll<nLEDs;ll+=1)
		if (NumType(CI[rr])==0) // if it is a cell used in clustering
			Kernel_SDCounter2[CurrentAngle_bin]+= 1
			Kernel_SDHists_summed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
			if (KernelList_SD[rr][ll]>SD_threshold) // ON kernels
				Kernel_SDHists[CurrentAngle_bin][CurrentIPLDepth_corrected][ll*2]+=1
				Kernel_SDCounter[CurrentAngle_bin][CurrentIPLDepth_corrected][ll*2,ll*2+1]+= 1
				Kernel_SDCounter_all[CurrentIPLDepth_corrected]+= 1	
				if (ll==0)
					Kernel_SDHists_Rsummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
				elseif (ll==1)
					Kernel_SDHists_Gsummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
				elseif (ll==2)
					Kernel_SDHists_Bsummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
				elseif (ll==3)
					Kernel_SDHists_Usummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1												
				endif
			elseif (KernelList_SD[rr][ll]<-SD_threshold) // OFF kernels			
				Kernel_SDHists[CurrentAngle_bin][CurrentIPLDepth_corrected][ll*2+1]-=-1
				Kernel_SDCounter[CurrentAngle_bin][CurrentIPLDepth_corrected][ll*2,ll*2+1]+= 1	
				Kernel_SDCounter_all[CurrentIPLDepth_corrected]+= 1		
				if (ll==0)
					Kernel_SDHists_Rsummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
				elseif (ll==1)
					Kernel_SDHists_Gsummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
				elseif (ll==2)
					Kernel_SDHists_Bsummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1
				elseif (ll==3)
					Kernel_SDHists_Usummed[CurrentAngle_bin][CurrentIPLDepth_corrected]+=1												
				endif
			endif
		endif
	endfor
endfor
Kernel_SDHists[][][]/=Kernel_SDCounter2[p]/nBins_depth
Kernel_SDHists_summed[][]/=Kernel_SDCounter2[p]/nBins_depth
Kernel_SDHists_Rsummed[][]/=Kernel_SDCounter2[p]/nBins_depth
Kernel_SDHists_Gsummed[][]/=Kernel_SDCounter2[p]/nBins_depth
Kernel_SDHists_Bsummed[][]/=Kernel_SDCounter2[p]/nBins_depth
Kernel_SDHists_Usummed[][]/=Kernel_SDCounter2[p]/nBins_depth

Wavestats/Q Kernel_SDCounter_all
//print "using", V_Sum, "/", nROIs*4, "entries (",100 * V_Sum/(nROIs*4),"%)"



if (AngleSmooth>0)
	make /o/n=(nBins_angle*3,nBins_depth,nLEDs*2) Kernel_SDHists_expanded = 0
	Kernel_SDHists_expanded[0,nBins_angle-1][][]=Kernel_SDHists[p][q][r]
	Kernel_SDHists_expanded[nBins_angle,2*nBins_angle-1][][]=Kernel_SDHists[p-nBins_angle][q][r]
	Kernel_SDHists_expanded[2*nBins_angle,3*nBins_angle-1][][]=Kernel_SDHists[p-nBins_angle*2][q][r]		
	Smooth /Dim=0 AngleSmooth, Kernel_SDHists_expanded
	Kernel_SDHists[][][]=Kernel_SDHists_expanded[p+nBins_angle][q][r]
	
	make /o/n=(nBins_angle*3,nBins_depth) Kernel_SDHists_expanded = 0
	Kernel_SDHists_expanded[0,nBins_angle-1][]=Kernel_SDHists_summed[p][q]
	Kernel_SDHists_expanded[nBins_angle,2*nBins_angle-1][]=Kernel_SDHists_summed[p-nBins_angle][q]
	Kernel_SDHists_expanded[2*nBins_angle,3*nBins_angle-1][]=Kernel_SDHists_summed[p-nBins_angle*2][q]		
	Smooth /Dim=0 AngleSmooth, Kernel_SDHists_expanded
	Kernel_SDHists_summed[][]=Kernel_SDHists_expanded[p+nBins_angle][q]
	Killwaves Kernel_SDHists_expanded
	
	make /o/n=(nBins_angle*3,nBins_depth) Kernel_SDHists_expanded = 0
	Kernel_SDHists_expanded[0,nBins_angle-1][]=Kernel_SDHists_Rsummed[p][q]
	Kernel_SDHists_expanded[nBins_angle,2*nBins_angle-1][]=Kernel_SDHists_Rsummed[p-nBins_angle][q]
	Kernel_SDHists_expanded[2*nBins_angle,3*nBins_angle-1][]=Kernel_SDHists_Rsummed[p-nBins_angle*2][q]		
	Smooth /Dim=0 AngleSmooth, Kernel_SDHists_expanded
	Kernel_SDHists_Rsummed[][]=Kernel_SDHists_expanded[p+nBins_angle][q]
	Killwaves Kernel_SDHists_expanded
	
	make /o/n=(nBins_angle*3,nBins_depth) Kernel_SDHists_expanded = 0
	Kernel_SDHists_expanded[0,nBins_angle-1][]=Kernel_SDHists_Gsummed[p][q]
	Kernel_SDHists_expanded[nBins_angle,2*nBins_angle-1][]=Kernel_SDHists_Gsummed[p-nBins_angle][q]
	Kernel_SDHists_expanded[2*nBins_angle,3*nBins_angle-1][]=Kernel_SDHists_Gsummed[p-nBins_angle*2][q]		
	Smooth /Dim=0 AngleSmooth, Kernel_SDHists_expanded
	Kernel_SDHists_Gsummed[][]=Kernel_SDHists_expanded[p+nBins_angle][q]
	Killwaves Kernel_SDHists_expanded
	
	make /o/n=(nBins_angle*3,nBins_depth) Kernel_SDHists_expanded = 0
	Kernel_SDHists_expanded[0,nBins_angle-1][]=Kernel_SDHists_Bsummed[p][q]
	Kernel_SDHists_expanded[nBins_angle,2*nBins_angle-1][]=Kernel_SDHists_Bsummed[p-nBins_angle][q]
	Kernel_SDHists_expanded[2*nBins_angle,3*nBins_angle-1][]=Kernel_SDHists_Bsummed[p-nBins_angle*2][q]		
	Smooth /Dim=0 AngleSmooth, Kernel_SDHists_expanded
	Kernel_SDHists_Bsummed[][]=Kernel_SDHists_expanded[p+nBins_angle][q]
	Killwaves Kernel_SDHists_expanded
	
	make /o/n=(nBins_angle*3,nBins_depth) Kernel_SDHists_expanded = 0
	Kernel_SDHists_expanded[0,nBins_angle-1][]=Kernel_SDHists_Usummed[p][q]
	Kernel_SDHists_expanded[nBins_angle,2*nBins_angle-1][]=Kernel_SDHists_Usummed[p-nBins_angle][q]
	Kernel_SDHists_expanded[2*nBins_angle,3*nBins_angle-1][]=Kernel_SDHists_Usummed[p-nBins_angle*2][q]		
	Smooth /Dim=0 AngleSmooth, Kernel_SDHists_expanded
	Kernel_SDHists_Usummed[][]=Kernel_SDHists_expanded[p+nBins_angle][q]
	Killwaves Kernel_SDHists_expanded
	
endif
if (DepthSmooth>0)
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists_summed	
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists_summed2
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists_Rsummed
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists_Gsummed
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists_Bsummed
	Smooth /Dim=1 DepthSmooth, Kernel_SDHists_Usummed			
endif
variable bb
for (bb=0;bb<nBins_angle;bb+=1)
	Kernel_SDHists_summed2[]+=Kernel_SDHists_summed[bb][p]/nBins_angle
	Kernel_SDHists_Rsummed2[]+=Kernel_SDHists_Rsummed[bb][p]/nBins_angle
	Kernel_SDHists_Gsummed2[]+=Kernel_SDHists_Gsummed[bb][p]/nBins_angle
	Kernel_SDHists_Bsummed2[]+=Kernel_SDHists_Bsummed[bb][p]/nBins_angle
	Kernel_SDHists_Usummed2[]+=Kernel_SDHists_Usummed[bb][p]/nBins_angle			
endfor

make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OnR = Kernel_SDHists[p][q][0]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OffR = Kernel_SDHists[p][q][1]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OnG = Kernel_SDHists[p][q][2]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OffG = Kernel_SDHists[p][q][3]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OnB = Kernel_SDHists[p][q][4]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OffB = Kernel_SDHists[p][q][5]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OnU = Kernel_SDHists[p][q][6]
make /o/n=(nBins_angle,nBins_depth) Kernel_SDHists_OffU = Kernel_SDHists[p][q][7]

make /o/n=(nBins_angle,nBins_depth) IPLMap_R_OnOff = 0
make /o/n=(nBins_angle,nBins_depth) IPLMap_G_OnOff = 0
make /o/n=(nBins_angle,nBins_depth) IPLMap_B_OnOff = 0
make /o/n=(nBins_angle,nBins_depth) IPLMap_U_OnOff = 0

IPLMap_R_OnOff[][]=((kernel_SDHists_OnR[p][q]) - (kernel_SDHists_OffR[p][q])) / ((kernel_SDHists_OnR[p][q]) + (kernel_SDHists_OffR[p][q]))
IPLMap_G_OnOff[][]=((kernel_SDHists_OnG[p][q]) - (kernel_SDHists_OffG[p][q])) / ((kernel_SDHists_OnG[p][q]) + (kernel_SDHists_OffG[p][q]))
IPLMap_B_OnOff[][]=((kernel_SDHists_OnB[p][q]) - (kernel_SDHists_OffB[p][q])) / ((kernel_SDHists_OnB[p][q]) + (kernel_SDHists_OffB[p][q]))
IPLMap_U_OnOff[][]=((kernel_SDHists_OnU[p][q]) - (kernel_SDHists_OffU[p][q])) / ((kernel_SDHists_OnU[p][q]) + (kernel_SDHists_OffU[p][q]))


Setscale x,-1,1,"pi radians" Kernel_SDHists
Setscale x,-1,1,"pi radians" Kernel_SDHists_OnR, Kernel_SDHists_OffR,IPLMap_R_OnOff
Setscale x,-1,1,"pi radians" Kernel_SDHists_OnG, Kernel_SDHists_OffG,IPLMap_G_OnOff
Setscale x,-1,1,"pi radians" Kernel_SDHists_OnB, Kernel_SDHists_OffB,IPLMap_B_OnOff
Setscale x,-1,1,"pi radians" Kernel_SDHists_OnU, Kernel_SDHists_OffU,IPLMap_U_OnOff


setscale y,0,100,"%" Kernel_SDHists, Kernel_SDHists_summed, Kernel_SDHists_summed2
setscale y,0,100,"%" Kernel_SDHists_Rsummed2, Kernel_SDHists_Gsummed2, Kernel_SDHists_Bsummed2, Kernel_SDHists_Usummed2
setscale x,0,100,"%" Kernel_SDCounter_all
setscale y,0,100,"%" Kernel_SDHists_OnR, Kernel_SDHists_OffR,IPLMap_R_OnOff
setscale y,0,100,"%" Kernel_SDHists_OnG, Kernel_SDHists_OffG,IPLMap_G_OnOff
setscale y,0,100,"%" Kernel_SDHists_OnB, Kernel_SDHists_OffB,IPLMap_B_OnOff
setscale y,0,100,"%" Kernel_SDHists_OnU, Kernel_SDHists_OffU,IPLMap_U_OnOff


make /o/n=(nBins_depth) currentFFT = 0
make /o/n=(nBins_angle,100-2*IPLE_edgecut) TempHist = Kernel_SDHists_Summed[p][q+IPLE_edgecut]

FFT/OUT=4/DEST=TempHist_FFT TempHist
currentFFT[]+=TempHist_FFT[0][p] 



// display

variable Display_range = 1

if (display_stuff==1)
	display /k=1 
	Appendimage /l=RedY /b=OnX Kernel_SDHists_OnR
	Appendimage /l=RedY /b=OffX Kernel_SDHists_OffR
	Appendimage /l=RedY /b=OnOffX IPLMap_R_OnOff	
	Appendimage /l=GreenY /b=OnX Kernel_SDHists_OnG
	Appendimage /l=GreenY /b=OffX Kernel_SDHists_OffG
	Appendimage /l=GreenY /b=OnOffX IPLMap_G_OnOff
	Appendimage /l=BlueY /b=OnX Kernel_SDHists_OnB
	Appendimage /l=BlueY /b=OffX Kernel_SDHists_OffB
	Appendimage /l=BlueY /b=OnOffX IPLMap_B_OnOff
	Appendimage /l=UVY /b=OnX Kernel_SDHists_OnU
	Appendimage /l=UVY /b=OffX Kernel_SDHists_OffU	
	Appendimage /l=UVY /b=OnOffX IPLMap_U_OnOff									
	
	ModifyImage Kernel_SDHists_OnR ctab= {0,Display_range,Red,0}
	ModifyImage Kernel_SDHists_OffR ctab= {0,Display_range,Red,0}
	ModifyImage IPLMap_R_OnOff ctab= {-1,1,BlueRedGreen,0}	
	ModifyImage Kernel_SDHists_OnG ctab= {0,Display_range,Green,0}
	ModifyImage Kernel_SDHists_OffG ctab= {0,Display_range,Green,0}
	ModifyImage IPLMap_G_OnOff ctab= {-1,1,BlueRedGreen,0}		
	ModifyImage Kernel_SDHists_OnB ctab= {0,Display_range,Blue,0}
	ModifyImage Kernel_SDHists_OffB ctab= {0,Display_range,Blue,0}
	ModifyImage IPLMap_B_OnOff ctab= {-1,1,BlueRedGreen,0}		
	ModifyImage Kernel_SDHists_OnU ctab= {0,Display_range,magenta,0}
	ModifyImage Kernel_SDHists_OffU ctab= {0,Display_range,magenta,0}	
	ModifyImage IPLMap_U_OnOff ctab= {-1,1,BlueRedGreen,0}	
	
	•ModifyGraph fSize=8,lblPos(OnX)=47,lblPos(OffX)=47,axisEnab(RedY)={0.8,1};DelayUpdate
	•ModifyGraph axisEnab(OnX)={0.05,0.5},axisEnab(OffX)={0.55,1};DelayUpdate
	•ModifyGraph axisEnab(GreenY)={0.55,0.75},axisEnab(BlueY)={0.3,0.5};DelayUpdate
	•ModifyGraph axisEnab(UVY)={0.05,0.25},freePos(RedY)={0,kwFraction};DelayUpdate
	•ModifyGraph freePos(OnX)={0,kwFraction},freePos(OffX)={0,kwFraction};DelayUpdate
	•ModifyGraph freePos(GreenY)={0,kwFraction},freePos(BlueY)={0,kwFraction};DelayUpdate
	•ModifyGraph freePos(UVY)={0,kwFraction};DelayUpdate
	•Label OnX "\\Z10 \U";DelayUpdate
	•Label OffX "\\Z10 \U"
	•Label OnOffX "\\Z10 \U";DelayUpdate
	
	•ModifyGraph noLabel(RedY)=1,noLabel(GreenY)=1,noLabel(BlueY)=1,noLabel(UVY)=1;DelayUpdate
	•ModifyGraph axThick(RedY)=0,axThick(GreenY)=0,axThick(BlueY)=0,axThick(UVY)=0;DelayUpdate
	•ModifyGraph lblPos=47
	ModifyGraph lblPos(GreenY)=0;DelayUpdate

	ModifyGraph noLabel(RedY)=2,noLabel(GreenY)=2,noLabel(BlueY)=2,noLabel(UVY)=2;DelayUpdate
	•ModifyGraph lblPos=47,axisEnab(OnX)={0,0.3},axisEnab(OffX)={0.35,0.65};DelayUpdate
	•ModifyGraph axisEnab(OnOffX)={0.7,1},freePos(OnOffX)={0,kwFraction};DelayUpdate
	•Label GreenY "";DelayUpdate
	•SetAxis OnX -1,1;DelayUpdate
	•SetAxis OffX -1,1;DelayUpdate
	•SetAxis OnOffX -1,1



	//
	
	Display /k=1
	String XONAxisname,XOFFAxisname,TraceName
	make /o/n=(nBins_depth) IPLPLotwave = -1
	setscale x,0,100,"%" IPLPlotwave

	for (bb=0;bb<nBins_angle;bb+=1)
		XONAxisName = "EyePosON"+Num2Str(bb)
		XOFFAxisName = "EyePosOFF"+Num2Str(bb)
		
		Appendtograph  /b=IPLY /l=$XONAxisName IPLPLotwave
		•ModifyGraph mode(IPLPLotwave#5)=7,useNegRGB(IPLPLotwave#5)=1;DelayUpdate
		•ModifyGraph hbFill(IPLPLotwave#5)=2,rgb(IPLPLotwave#5)=(52224,52224,52224);DelayUpdate
		•ModifyGraph negRGB(IPLPLotwave#5)=(52224,52224,52224)
		
		
		Tracename = "Kernel_SDHists_OnR#"+Num2Str(bb*8)
		if (bb==0)
			Tracename = "Kernel_SDHists_OnR"
		endif
		Appendtograph /b=IPLY /l=$XONAxisName Kernel_SDHists_OnR[bb][]
		Tracename = "Kernel_SDHists_OffR#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OffR"
		endif
		Appendtograph /b=IPLY /l=$XOFFAxisName Kernel_SDHists_OffR[bb][]
		Tracename = "Kernel_SDHists_OnG#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OnG"
		endif
		Appendtograph /b=IPLY /l=$XONAxisName Kernel_SDHists_OnG[bb][]
		ModifyGraph rgb($Tracename)=(0,2^16-1,0), offset($tracename)={-offsetScale*1,0}
		Tracename = "Kernel_SDHists_OffG#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OffG"
		endif
		Appendtograph /b=IPLY /l=$XOFFAxisName Kernel_SDHists_OffG[bb][]
		ModifyGraph rgb($Tracename)=(0,2^16-1,0), offset($tracename)={-offsetScale*1,0}
		Tracename = "Kernel_SDHists_OnB#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OnB"
		endif
		Appendtograph /b=IPLY /l=$XONAxisName Kernel_SDHists_OnB[bb][]
		ModifyGraph rgb($Tracename)=(0,0,2^16-1), offset($tracename)={-offsetScale*2,0}
		Tracename = "Kernel_SDHists_OffB#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OffB"
		endif
		Appendtograph /b=IPLY /l=$XOFFAxisName Kernel_SDHists_OffB[bb][]
		ModifyGraph rgb($Tracename)=(0,0,2^16-1), offset($tracename)={-offsetScale*2,0}
		Tracename = "Kernel_SDHists_OnU#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OnU"
		endif
		Appendtograph /b=IPLY /l=$XONAxisName Kernel_SDHists_OnU[bb][]
		ModifyGraph rgb($Tracename)=(29440,0,58880), offset($tracename)={-offsetScale*3,0}	
		Tracename = "Kernel_SDHists_OffU#"+Num2Str(bb)
		if (bb==0)
			Tracename = "Kernel_SDHists_OffU"
		endif
		Appendtograph /b=IPLY /l=$XOFFAxisName Kernel_SDHists_OffU[bb][]
		ModifyGraph rgb($Tracename)=(29440,0,58880), offset($tracename)={-offsetScale*3,0}
	
		variable XPlotFrom = 0.05 + bb * 0.95/nBins_angle
		variable XPLotTo = XPlotFrom + 0.8/nBins_angle
		
		
		ModifyGraph axisEnab($XONAxisName)={XPLotfrom,XPlotTo}
		ModifyGraph axisEnab($XOFFAxisName)={XPLotfrom,XPlotTo}
		SetAxis $XONAxisName -1,1
		SetAxis $XOFFAxisName 1,-1
		ModifyGraph zero($XONAxisName)=2
		ModifyGraph noLabel($XOFFAxisName)=2,axThick($XOFFAxisName)=0
		
	
	endfor
	
	•ModifyGraph mode=0,marker=16;DelayUpdate
	•ModifyGraph msize=1.5
	ModifyGraph lsize=1.5
	
	for (bb=0;bb<nBins_angle;bb+=1)
		Tracename = "IPLPLotwave#"+Num2Str(bb)
		if (bb==0)
			Tracename = "IPLPLotwave#"
		endif
		
		•ModifyGraph mode($Tracename)=7,useNegRGB($Tracename)=1;DelayUpdate
		•ModifyGraph hbFill($Tracename)=2,rgb($Tracename)=(52224,52224,52224);DelayUpdate
		•ModifyGraph negRGB($Tracename)=(52224,52224,52224)
	
	endfor
	
	ModifyGraph swapXY=1

	ModifyGraph fSize=8,freePos={0,kwFraction}

	•ModifyGraph lblPos(IPLY)=47,axisEnab(IPLY)={0.05,1};DelayUpdate
	•Label IPLY "\\Z10IPL depth (\U)"
	SetAxis IPLY 0,100



endif



end