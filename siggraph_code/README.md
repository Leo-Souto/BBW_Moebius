Matlab source code for SGP 2020 Submission ID 1026
Handling conflicts on the sphere with Moebius transformations

System Requirements
- MATLAB 2014a or more recent version
- Image Processing Toolbox
- Optimization Toolbox

Installation instructions
- Clone gptoolbox repository at: https://github.com/alecjacobson/gptoolbox
- Add this folder and subfolders to MATLAB path
- Add gptoolbox folder and subfolders to MATLAB path

Using the User interface
- On Matlab's command window, issue the following command:
run biharmonic_interface.m
- Load the image and add the desired handles by left-clicking over the image. Note: some users have experienced misplacement of handles when they click on the image. This problem is solved by maximizing the window.
- To finish a closed cage handle right-click the image
- Press next 
- Rectify all lines by clicking in rectify lines or
- Rectify current line by clicking in rectify current line or
- Select control point by clicking over it (including point handle arrows)
- Change control point position by clicking over the image
- Undo the selected control point transformation by clicking in the respective button
- Click save image to generate the final result

This code has been tested in the following settings:
- Linux (Ubuntu 15.10) with Matlab 2014a
- Mac OS (High Sierra) with Matlab 2017a
