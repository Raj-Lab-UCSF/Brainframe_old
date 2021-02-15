function input_struct = brainframe_inputs_mouse(matpath,varargin)

%Setting the path to load defaults; set matpath to [] if using default
if mod(nargin,2) == 0 || isempty(matpath)
%     matpath = [cd filesep 'Data'];
    matpath = cd;
end

%Loading the data necessary to create this default input_struct:
input_struct = struct;
load([matpath filesep 'brainframe_defaultMouse_datinput.mat'],...
    'brainat','conmat','pathology');


%Setting default parameter values
%Specific fields are filled one at a time in the below script
%Basic fields required for any element of brainframe to run:

%Voxel or region binary flag, enter 0 or 1 respectively
voxUreg_ = 1;

%Setting 3D brain atlas with regional voxel IDs, mouse default is slightly modified AIBS CCF
brain_atlas_ = brainat;

%Setting image background color, only takes 'k','w','other', where other is light gray
bgcolor_ = 'k';

%Binary flag for saving & closing image with on axis views [1], or opening GUI for image manipulation [0]
savenclose_ = 0;

%Field to set image file name if auto-saving
img_labels_ = 'yourfilename';

%Field to set image file format
img_format_ = 'png';

%Field for entering your data to be visualized over the space of the brain
%This field is a vector of values for per-region data that is n-regions long
%This field is either a 3D matrix the same size as or a 1D vector with the same number of elements as the brain atlas
%Per region semi-quantitative tau pathology data is the default example form mice from a figure in Iba, et al., 2013 (see bibliography in help file)
pathology(isnan(pathology)) = 0;
data_ = pathology(:,1);

%Relevant For Per Voxel Visualizations Only:

%Number of evenly spaced bins for heatmap visualization of per voxel data
nbin_ = 1;

%Relevant For Per Region Visualizations: 

%Setting number of regions
nreg_ = sum(unique(brain_atlas_)>0);

%Setting regions groups, here all the same group
region_groups_ = ones(sum(unique(brain_atlas_)>0),1);

%Binary flag for whether to draw a sphere [1] or point clouds [0] for per-region visualizations
sphere_ = 0;

%Setting number of points in per sphere visualization
sphere_npts_ = 35;

%Index 1: binary flag for centered [1] or diffuse [0] point clouds
%Index 2: Degree of centering, only turned on if [1] is in index 1
centered_ = [1 2];


%Manipulable Fields Relevant For Both Per-Region & Per Voxel Visualizations:

%Field for setting colormap, based on nbins for per-voxel & region_groups for per-region visualizations
cmap_ = hsv(length(unique(region_groups_)));

%Multiplier for size & density of all point clouds or radius of all spheres
xfac_ = 1;

%Point size specification (uses ptcloud(), scatter3() behaves similarly)
pointsize_ = 1;


%Fields Relevant For Connectivity Visualizations:

%Binary flag for visualizing interregional connectivity data, only visually works with per-region visualizations
iscon_ = 0;

%Field for setting connectivity matrix, must have same number of regions as brain_atlas
conmat_ = conmat;

%Multiplier for number of fibers from conmat
con_rescale_ = 0.01;

%Fiber width setting (uses plot3() for visualizing)
con_width_ = 0.01;

%Field for setting connections into regions groups based on ROI of origin, functions similarly to region_groups above
con_regiongroups_ = ones(sum(unique(brain_atlas_)>0),1);

%Field for setting connectivity colormap
con_cmap_ = lines(length(unique(region_groups_)));

%Field for setting degree of curvature of ellipse connectivity visualizations
con_arch_ = 0.5;

%Setting the width and length of arrows on fibers indicating direction (pre to post synpatic ROI)
conarrow_WL_ = [0.5 0.8];

%Creating input parser
ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validNonnegative = @(x) isnumeric(x) && all(x(:)>=0);
validBoolean = @(x) isscalar(x) && (x==0 || x==1);
validChar = @(x) ischar(x);
validImg = @(x) strcmp(x,'png') || strcmp(x,'jpg') || strcmp(x,'tiffn');
validCentered = @(x) validBoolean(x(1)) && validScalar(x(2)) && length(x) == 2;
validConarrow = @(x) validNonnegative(x) && length(x) == 2;

addParameter(ip, 'voxUreg', voxUreg_, validBoolean);
addParameter(ip, 'brain_atlas', brain_atlas_, validNonnegative);
addParameter(ip, 'bgcolor', bgcolor_, validChar);
addParameter(ip, 'savenclose', savenclose_, validBoolean);
addParameter(ip, 'img_labels', img_labels_, validChar);
addParameter(ip, 'img_format', img_format_, validImg);
addParameter(ip, 'data', data_, validNonnegative);
addParameter(ip, 'nbin', nbin_, validScalar);
addParameter(ip, 'nreg', nreg_, validScalar);
addParameter(ip, 'region_groups', region_groups_, validNonnegative);
addParameter(ip, 'sphere', sphere_, validBoolean);
addParameter(ip, 'sphere_npts', sphere_npts_, validScalar);
addParameter(ip, 'centered', centered_, validCentered);
addParameter(ip, 'cmap', cmap_, validNonnegative);
addParameter(ip, 'xfac', xfac_, validScalar);
addParameter(ip, 'pointsize', pointsize_, validScalar);
addParameter(ip, 'iscon', iscon_, validBoolean);
addParameter(ip, 'conmat', conmat_, validNonnegative);
addParameter(ip, 'con_rescale', con_rescale_, validScalar);
addParameter(ip, 'con_width', con_width_, validScalar);
addParameter(ip, 'con_regiongroups', con_regiongroups_, validNonnegative);
addParameter(ip, 'con_cmap', con_cmap_, validNonnegative);
addParameter(ip, 'con_arch', con_arch_, validScalar);
addParameter(ip, 'conarrow_WL', conarrow_WL_, validConarrow);
parse(ip, varargin{:});

input_struct.voxUreg = ip.Results.voxUreg;
input_struct.brain_atlas = ip.Results.brain_atlas;
input_struct.bgcolor = ip.Results.bgcolor;
input_struct.savenclose = ip.Results.savenclose;
input_struct.img_labels = ip.Results.img_labels;
input_struct.img_format = ip.Results.img_format;
input_struct.data = ip.Results.data;
input_struct.nbin = ip.Results.nbin;
input_struct.nreg = ip.Results.nreg;
input_struct.region_groups = ip.Results.region_groups;
input_struct.sphere = ip.Results.sphere;
input_struct.sphere_npts = ip.Results.sphere_npts;
input_struct.centered = ip.Results.centered;
input_struct.cmap = ip.Results.cmap;
input_struct.xfac = ip.Results.xfac;
input_struct.pointsize = ip.Results.pointsize;
input_struct.iscon = ip.Results.iscon;
input_struct.conmat = ip.Results.conmat;
input_struct.con_rescale = ip.Results.con_rescale;
input_struct.con_width = ip.Results.con_width;
input_struct.con_regiongroups = ip.Results.con_regiongroups;
input_struct.con_cmap = ip.Results.con_cmap;
input_struct.con_arch = ip.Results.con_arch;
input_struct.conarrow_WL = ip.Results.conarrow_WL;

end


