Matlab source code for SIGGRAPH Submission ID 432
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
run mesh_interface.m
- Choose the quality of the mesh by changing the quality slider
- Press next
- Load the image and add the desired handles by left-clicking over the image
- To finish a closed cage handle right-click the image
- To add a conformal region, click in the image and drag the cross over the region
- Press next to compute the BBW
- Check if weights are good as desired
- Press next
- Select control point by clicking over it (including point handle arrows)
- Change control point position by clicking over the image
- Undo the selected control point transformation by clicking in the respective button
- Rectify all lines by clicking in rectify lines
- Rectify current line by clicking in rectify current line
- Click save image to generate the final result


Some known bugs that will be fixed in the future
- Sometimes the mesh created lack points on the poles, change slightly the quality slider to correct
- The delete handle button is not working properly with closed cage handles
- Sometimes the middle point of bone handles disappear after rectification


This code has been tested in the following settings:
- Linux (Ubuntu 15.10) with Matlab 2014a
- Mac OS (High Sierra) with Matlab 2017a
