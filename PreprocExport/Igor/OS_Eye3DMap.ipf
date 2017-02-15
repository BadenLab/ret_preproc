#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function OS_Eye3DMap()

wave wParamsNum
wave positions
wave ROIs

variable Xpos =  FindDimLabel(wParamsNum,0,"XCoord_um")
variable Ypos =  FindDimLabel(wParamsNum,0,"YCoord_um")
variable Zpos =  FindDimLabel(wParamsNum,0,"ZCoord_um")

variable Zoom =  FindDimLabel(wParamsNum,0,"Zoom")
variable Angle =  FindDimLabel(wParamsNum,0,"Angle_deg")
variable X_UserOffset =  FindDimLabel(wParamsNum,0,"User_XOffset_V")
variable Y_UserOffset =  FindDimLabel(wParamsNum,0,"User_YOffset_V")

variable nLines = wParamsNum[18]
variable xPixelsInd = FindDimLabel(wParamsNum,0,"User_dxPix" )
variable yPixelsInd = FindDimLabel(wParamsNum,0,"User_dyPix" )
variable realPixDurInd = FindDimLabel(wParamsNum,0,"RealPixDur" )
variable lineDur = (wParamsNum[xPixelsInd] *  wParamsNum[realPixDurInd]) * 10^-6

variable nROIs = Dimsize(positions,0)

//////////////////////////////////////







end