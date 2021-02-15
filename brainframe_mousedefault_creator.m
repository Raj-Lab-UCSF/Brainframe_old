%Loading the data necessary to create this default input_struct:

matpath = cd;
load([matpath filesep 'brainframe_defaultMouse_datinput.mat']);


%Creating the default mouse input struct
%Specific fields are filled one at a time in the below script
%Basic fields required for any element of brainframe to run:

%Voxel or region binary flag, enter 0 or 1 respectively
input_struct.voxUreg = 1;

%Setting 3D brain atlas with regional voxel IDs, mouse default is slightly modified AIBS CCF
input_struct.brain_atlas = brainat;

%Setting image background color, only takes 'k','w','other', where other is light gray
input_struct.bgcolor = 'k';

%Binary flag for saving & closing image with on axis views [1], or opening GUI for image manipulation [0]
input_struct.savenclose = 0;

%Field to set image file name if auto-saving
input_struct.img_labels = 'yourfilename';

%Field to set image file format
input_struct.img_format = 'png';

%Field for entering your data to be visualized over the space of the brain
%This field is a vector of values for per-region data that is n-regions long
%This field is either a 3D matrix the same size as or a 1D vector with the same number of elements as the brain atlas
%Per region semi-quantitative tau pathology data is the default example form mice from a figure in Iba, et al., 2013 (see bibliography in help file)
pathology(isnan(pathology)) = 0;
input_struct.data = pathology(:,1);

%Relevant For Per Voxel Visualizations Only:

%Number of evenly spaced bins for heatmap visualization of per voxel data
input_struct.nbin = 1;


%Relevant For Per Region Visualizations: 

%Setting number of regions
input_struct.nreg = sum(unique(input_struct.brain_atlas)>0);

%Setting regions groups, here all the same group
input_struct.region_groups = ones(sum(unique(input_struct.brain_atlas)>0),1);

%Binary flag for whether to draw a sphere [1] or point clouds [0] for per-region visualizations
input_struct.sphere = 0;

%Setting number of points in per sphere visualization
input_struct.sphere_npts = 35;

%Index 1: binary flag for centered [1] or diffuse [0] point clouds
%Index 2: Degree of centering, only turned on if [1] is in index 1
input_struct.centered = [1 2];


%Manipulable Fields Relevant For Both Per-Region & Per Voxel Visualizations:

%Field for setting colormap, based on nbins for per-voxel & region_groups for per-region visualizations
input_struct.cmap = hsv(length(unique(input_struct.region_groups)));

%Multiplier for size & density of all point clouds or radius of all spheres
input_struct.xfac = 1;

%Point size specification (uses ptcloud(), scatter3() behaves similarly)
input_struct.pointsize = 1;


%Fields Relevant For Connectivity Visualizations:

%Binary flag for visualizing interregional connectivity data, only visually works with per-region visualizations
input_struct.iscon = 0;

%Field for setting connectivity matrix, must have same number of regions as brain_atlas
input_struct.conmat = conmat;

%Multiplier for number of fibers from conmat
input_struct.con_rescale = 0.01;

%Fiber width setting (uses plot3() for visualizing)
input_struct.con_width = 0.01;

%Field for setting connections into regions groups based on ROI of origin, functions similarly to region_groups above
input_struct.con_regiongroups = ones(sum(unique(input_struct.brain_atlas)>0),1);

%Field for setting connectivity colormap
input_struct.con_cmap = lines(length(unique(input_struct.region_groups)));

%Field for setting degree of curvature of ellipse connectivity visualizations
input_struct.con_arch = 0.5;

%Setting the width and length of arrows on fibers indicating direction (pre to post synpatic ROI)
input_struct.conarrow_WL = [0.5 0.8];


%Saving this default input_struct:

save([matpath filesep 'default_mouse.mat'],'input_struct');

