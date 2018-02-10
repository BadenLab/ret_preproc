#pragma rtGlobals=3		// Use modern global access method and strict wave access.

#include "OS_ParameterTable"
#include "OS_DetrendStack"
#include "OS_ManualROI"
#include "OS_AutoRoiByCorr"
#include "OS_TracesAndTriggers"
#include "OS_BasicAveraging"
#include "OS_hdf5Export"
#include "OS_LaunchCellLab"
#include "OS_STRFs"
#include "OS_EventFinder"
#include "OS_hdf5Import"
#include "OS_LineScanFormat"
#include "OS_LED_Noise"
#include "OS_Clustering"
#include "OS_Skittles"
#include "OS_KernelfromROI" 
#include "OS_STRFs" // NOT ADDED YET


//----------------------------------------------------------------------------------------------------------------------
Menu "ScanM", dynamic
	"-"
	" Open OS GUI",	/Q, 	OS_GUI()
	"-"	
End
//----------------------------------------------------------------------------------------------------------------------


function OS_GUI()
	NewPanel /N=OfficialScripts /k=1 /W=(200,100,450,660)
	ShowTools/A
	SetDrawLayer UserBack

	SetDrawEnv fstyle= 1
	DrawText 24,36,"(Step 0: Linescan only)"
	SetDrawEnv fstyle= 1
	DrawText 24,36+54,"Step 1: Parameter Table"
	SetDrawEnv fstyle= 1
	DrawText 24,90+54,"Step 2: Pre-formatting"
	SetDrawEnv fstyle= 1
	DrawText 24,149+54,"Step 3: ROI placement"
	SetDrawEnv fstyle= 1
	DrawText 24,272+54,"Step 4: Extract Traces and Triggers"
	SetDrawEnv fstyle= 1
	DrawText 24,334+54,"Step 5a: Further optional processes"
	SetDrawEnv fstyle= 1	
	DrawText 24,454+54,"Step 6: Database Export/Import (hdf5)"
	Button step0,pos={99,39},size={100,26},proc=OS_GUI_Buttonpress,title="Process linescan"
	Button step1a,pos={78,39+54},size={107,26},proc=OS_GUI_Buttonpress,title="Make / Show"
	Button step1b,pos={192,39+54},size={34,26},proc=OS_GUI_Buttonpress,title="Kill"	
	Button step2a,pos={78,94+54},size={71,26},proc=OS_GUI_Buttonpress,title="Standard"
	Button step2b,pos={154,94+54},size={71,26},proc=OS_GUI_Buttonpress,title="Ratiom."
	Button step3a1,pos={78,155+54},size={43,20},proc=OS_GUI_Buttonpress,title="Manual"
	Button step3a2,pos={130,155+54},size={43,20},proc=OS_GUI_Buttonpress,title="Apply"
	Button step3a3,pos={181,155+54},size={43,20},proc=OS_GUI_Buttonpress,title="Pixels"
	Button step3a4,pos={78,179+54},size={147,20},proc=OS_GUI_Buttonpress,title="Use existing SARFIA Mask"	
	Button step3b,pos={78,203+54},size={147,20},proc=OS_GUI_Buttonpress,title="Autom. by Correlation"
	Button step3c,pos={78,228+54},size={147,20},proc=OS_GUI_Buttonpress,title="Autom. CellLab"
	Button step4,pos={78,278+54},size={147,26},proc=OS_GUI_Buttonpress,title="Traces and Triggers"
	Button step5a,pos={78,341+54},size={43,26},proc=OS_GUI_Buttonpress,title="Ave"
	Button step5b,pos={130,341+54},size={43,26},proc=OS_GUI_Buttonpress,title="Events"			
	Button step5c,pos={181,341+54},size={43,26},proc=OS_GUI_Buttonpress,title="Kernels"	
	Button step5d,pos={78,371+54},size={43,26},proc=OS_GUI_Buttonpress,title=" Cluster "			
	Button step5e,pos={130,371+54},size={43,26},proc=OS_GUI_Buttonpress,title=" ROI-K"	
	Button step5f,pos={181,371+54},size={43,26},proc=OS_GUI_Buttonpress,title=" K-Map "	

	Button step5g,pos={78,401+54},size={71,26},proc=OS_GUI_Buttonpress,title=" Sweep "			
	Button step5h,pos={154,401+54},size={71,26},proc=OS_GUI_Buttonpress,title=" Swoosh "
	
	Button step6a,pos={78,462+54},size={71,26},proc=OS_GUI_Buttonpress,title="Export"
	Button step6b,pos={154,462+54},size={71,26},proc=OS_GUI_Buttonpress,title="Import"	
	
	HideTools/A
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Function OS_GUI_Buttonpress(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			strswitch (ba.ctrlName)
				case "step0":
					OS_LineScanFormat()
					break
				case "step1a":
					OS_ParameterTable()
					break
				case "step1b":
					OS_ParameterTable_Kill()
					break					
				case "step2a":
					OS_DetrendStack()
					break
				case "step2b":
					OS_DetrendRatiometric()
					break					
				case "step3a1":
					OS_CallManualROI()
					break
				case "step3a2":
					OS_ApplyManualRoi()
					break	
				case "step3a3":
					OS_monoPixelApply()
					break						
				case "step3a4":
					OS_CloneSarfiaRoi()
					break																		
				case "step3b":
					OS_AutoRoiByCorr()
					break
				case "step3c":
					OS_LaunchCellLab()
					break
				case "step4":
					OS_TracesAndTriggers()
					break					
				case "step5a":
					OS_BasicAveraging()
					break
				case "step5b":
					OS_EventFinder()
					break					
				case "step5c":
					OS_LED_Noise()//OS_STRFs() // replaced by kernel script for 4LED Noise
					break
				case "step5d":
					OS_Clustering()
					break
				case "step5e":
					OS_KernelfromROI()
					break
				case "step5f":
					OS_IPLKernels()
					break	
				case "step5g":
					OS_SkittlesSweep()
					break		
				case "step5h":
					OS_SkittlesSwooshMap()
					break																					
				case "step6a":
					OS_hdf5Export()
					break										
				case "step6b":
					OS_hdf5Import("")
					break
			endswitch
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
