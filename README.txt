Program name: 
	An Enhanced Spatiotemporal Data Fusion Method Integrating Stepwise Conversion Coefficient Calculation and Flexible Pixel Prediction Strategy (E2STARFM)
	Version 1:
	School of Geography and Information Engineering, China University of Geosciences (Wuhan)
	August 21, 2024
	

Purpose: 
	The purpose of this method is to fuse two pairs of fine-resolution and coarse-resolution images from different sensors at Tm and Tn, and a coarse-resolution image at Tp, to produce a fine-resolution image at Tp. 
	

Usage: 
	There are several parameters for performance tuning:
		w : The haif window size, if 25, the window size is 25*2+1=51 fine pixels.
    		DN_min & DN_max : The range of DN value of the image, If byte, 0 and 255.
    		background_value : The value of background and missng pixels in both MODIS and Landsat images.
    		d_param : The free parameter. Parameter d may be slightly different for different sensors.
    		patch_long : The size of each block, if process whole ETM scene, set 500-1000.
    		pre_file : The path for storing temporary files.


Inputs: 
		FileName1 : fine-resolution image at Tm.
		FileName2 : coarse-resolution image at Tm.
   		FileName3 : fine-resolution image at Tn.
    		FileName4 : coarse-resolution image at Tn.
    		FileName5 : coarse-resolution image at Tp.


Output: 
	fine-resolution image at Tp.


Contact：
	If you have any question, please contact Chuchu She (cshe743248@cug.edu.cn) or submit a issue.
