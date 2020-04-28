*** Note to self: KEEP POPULATING THIS! ***

Repository for Axial and Pinon Flat SCTA projects
________________________________________________________________________________
> process_Axial_ver2.m reads Axial SCTA data from Google Drive (preferred) or
from IRIS and extracts acceleration values for calibrations.

Outputs:
    ~/Documents/SCTA/calibrations/Axial/axialdata.mat
        dataDec1 - t,x,y,z,a,T data at 1 Hz
        dataDec100 - t,x,y,z,a,T data at 1/100 Hz
        flipInfoAll - various data from calibrations

Plots:
    ~/Documents/SCTA/calibrations/Axial/process_Axial_ver2.tiff
        scatter plot of calibrations vs. time
________________________________________________________________________________
> process_PinonFlat_ver2.m reads Pinon Flat SCTA data from Google Drive (preferred)
and extracts acceleration values for calibrations.

Outputs:
    ~/Documents/PinonFlat/calibrations/PF/PFdata.mat
        dataDec1 - t,x,y,z,a,T data at 1 Hz
        dataDec100 - t,x,y,z,a,T data at 1/100 Hz
        flipInfoAll - various data from calibrations

Plots:
    ~/Documents/PinonFlat/calibrations/PF/process_PinonFlat_ver2.tiff
        scatter plot of calibrations vs. time
________________________________________________________________________________
> longterm_tiltplots.m generates continuous time series of SCTA data at 1 sample/minute
and 1 sample/hour

Inputs:
    ~/Documents/SCTA/calibrations/Axial/axialdata.mat

Outputs:
    ~/Documents/SCTA/calibrations/Axial/axialdata_hr.mat
        data_hr - t,x,y,z,T data at 1 sample/hour
    ~/Documents/SCTA/calibrations/Axial/axialdata_min.mat
        data_min - t,x,y,z,T data at 1 sample/minute
________________________________________________________________________________
> stitchAxial_8Hz_ver2.m removes calibration intervals from continuous tilt data,
corrects for offsets, and applies linear drift corrections based on calibrations

Inputs:
    ~/Documents/SCTA/calibration/Axial/axialdata_hr.mat
    ~/Documents/SCTA/calibration/Axial/axialdata_min.mat

Outputs
    ~/Documents/SCTA/calibration/Axial/axialstitch_hr.mat
        stitch_hr - t,x,y,z,T data at 1 sample/hour
    ~/Documents/SCTA/calibration/Axial/axialstitch_min.mat
        stitch_min - t,x,y,z,T data at 1 sample/minute
________________________________________________________________________________
> drift_spandrift.m quantifies the drift of the x and y SCTA channels, as well
as the drift of the span of those channels

Inputs:
    ~/Documents/SCTA/calibrations/Axial/axialdata.mat
        flipInfoAll - various data from calibrations

Plots:
    ~/Documents/SCTA/calibrations/Axial/span/xspandrift.tiff
        time series of + and - x calibrations and span (normalized to fit on axes)
        with linear best fits and slope values printed
    ~/Documents/SCTA/calibrations/Axial/span/yspandrift.tiff
        time series of + and - y calibrations and span (normalized to fit on axes)
        with linear best fits and slope values printed
