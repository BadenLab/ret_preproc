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


//----------------------------------------------------------------------------------------------------------------------
Menu "ScanM", dynamic
	"-"
	" Open OS GUI",	/Q, 	OS_GUI()
	"-"	
End
//----------------------------------------------------------------------------------------------------------------------


function OS_GUI()
	NewPanel /N=OfficialScripts /k=1 /W=(200,100,450,650)
	ShowTools/A
	SetDrawLayer UserBack

	SetDrawEnv fstyle= 1
	DrawText 24,36,"(Step 0: Linescan only)"
	SetDrawEnv fstyle= 1
	DrawText 24,36+54,"Step 1: Generate Parameter Table"
	SetDrawEnv fstyle= 1
	DrawText 24,90+54,"Step 2: Detrending"
	SetDrawEnv fstyle= 1
	DrawText 24,149+54,"Step 3: ROI placement"
	SetDrawEnv fstyle= 1
	DrawText 24,272+54,"Step 4: Extract Traces and Triggers"
	SetDrawEnv fstyle= 1
	DrawText 24,334+54,"Step 5a: Further optional processes"
	SetDrawEnv fstyle= 1	
	DrawText 24,424+54,"Step 6: Database Export/Import (hdf5)"
	Button step0,pos={99,39},size={100,26},proc=OS_GUI_Buttonpress,title="Process linescan"
	Button step1,pos={78,39+54},size={147,26},proc=OS_GUI_Buttonpress,title="Make/Show Parameters"
	Button step2a,pos={78,94+54},size={71,26},proc=OS_GUI_Buttonpress,title="One Channel"
	Button step2b,pos={154,94+54},size={71,26},proc=OS_GUI_Buttonpress,title="Ratiometric"
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
	Button step5d,pos={78,371+54},size={71,26},proc=OS_GUI_Buttonpress,title=" QuickCluster "			
	Button step5e,pos={154,371+54},size={71,26},proc=OS_GUI_Buttonpress,title=" KernelMap "	
	Button step6a,pos={78,432+54},size={71,26},proc=OS_GUI_Buttonpress,title="Export"
	Button step6b,pos={154,432+54},size={71,26},proc=OS_GUI_Buttonpress,title="Import"	
	
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
				case "step1":
					OS_ParameterTable()
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
					OS_IPLKernels()
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
