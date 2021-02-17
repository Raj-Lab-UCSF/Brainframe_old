# Brainframe
[Chris Mezias](https://github.com/chm2062) & [Justin Torok](https://github.com/justin-torok)

Matlab functionality for creating 3D human & mouse brain renderings, with many options. This is a tool for visualizing per-voxel or per-region data, as well as connectome data, in a 3D surface representation of a brain. This function has been tested for and works using mice and human brain atlases, but may well work with other species' brains. Further help including use case examples can be found in the script `brainframe_Help.m`.

To begin using the Brainframe package, please clone or download the entire repository into a folder accessible to Matlab on your computer. Rearranging folder contents before getting comfortable with use could lead to package functionality breaking.

## 1. Functions, Scripts, & Data Files
Below is a short description of the various functions, scripts, and data files included in the Brainframe package.

### Functions
- `brainframe.m`: the rendering function, which takes an input_struct, with fields described in the below **Fields of input_struct** subsection.

- `brainframe_inputs_human`: creates a human input_struct for `brainframe.m` to render. There are two types of input this functions takes. The first is a filepath which the function loads from and saves the input_struct into. The second set of inputs follows the `varargin` style of Matlab. Only specifying the filepath will result in the default input_struct being created. Please the below code snippet for an example:

	> matpath = cd; %Change this to alter the path you load from
	> input_struct = brainframe_inputs_human(matpath);

	- The above will generate a human input_struct which can then be fed to `brainframe.m` as follows, generating the image seen below:

	> brainframe(input_struct)

	![human_default](/ExampleImages/brainframe_Help_01.png)

- `brainframe_inputs_mouse`: creates a mouse input_struct. This function operates exactly the same way as the human function, with exactly the same input_struct fields. Please see below for an example:

	> matpath = cd; %Change this to alter the path you load from
	> input_struct = brainframe_inputs_mouse(matpath);
	> brainframe(input_struct)
	> view([-1 0 0]);

	![mouse_default](/ExampleImages/brainframe_Help_02.png)

- `arrow3`: a widely used Mathworks exchange function, which the authors of Brainframe did not create. This function is only called internally by `brainframe.m` during rendering and it is not necessary for the user to interact with it. Please refer to the following link for more information about this function: [arrow3.m](https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3)

### Scripts
- `brainframe_Help`: an extensively commented example script which goes through the basic syntax of using the Brainframe package. This script also contains almost a dozen use case examples, and if run will generate images following the code examples.

- `brainframe_humandefault_creator.m`: this script is an extensively commented example demonstrating both the default human input_struct field values and how to change or declare input_struct fields directly by using an input_struct.fieldname syntax. This can be found in the **ExampleDefaults** folder.

- `brainframe_mousedefault_creator.m`: this script is an extensively commented example demonstrating both the default mouse input_struct field values and how to change or declare input_struct fields directly by using an input_struct.fieldname syntax. This can be found in the **ExampleDefaults** folder.

### Data Files
- **brainframe_defaultHuman_datinput.mat**: this function contains the volumetric 86-region Desikan atlas, 86-region Desikan connectome, and an 86 element vector with example tau pathology data, per-region, from an average of human patients' tau PET data from ADNI. It is necessary to have this data file for both the `brainframe_humandefault_creator.m` script and `brainframe_inputs_human` function to create the standard human default input_struct.

- **brainframe_defaultMouse_datinput.mat**: this function contains the volumetric 426-region AIBS atlas, 426-region AIBS connectome, and a 426 element vector with example tau pathology data, per-region, from an average of mice 1 month after tau injection. It is necessary to have this data file for both the `brainframe_mousedefault_creator.m` script and `brainframe_inputs_mouse` function to create the standard human default input_struct.

- **default_human.mat**: this data file is the default input_struct created by both the `brainframe_humandefault_creator.m` script and `brainframe_inputs_human` function when no Name and Value pair arguments are specified.

- **default_mouse.mat**: this data file is the default input_struct created by both the `brainframe_mousedefault_creator.m` script and `brainframe_inputs_mouse` function when no Name and Value pair arguments are specified.

- **PerVox_ExampleData.mat**: this is the example per-voxel data that can be rendered. Please see the **Brief Examples** subsection below and lines 207-216 in the script `brainframe_Help.m` for further information and example renderings.

## 2. Fields of input_struct
This section describes the fields of input_struct, the structure object which `brainframe.m` takes for brain rendering. Modifying the fields of input_struct will change the rendering. You can change these fields directly in an already specified input_struct. Alternately, you can create a new input_struct with either `brainframe_inputs_human` or `brainframe_inputs_mouse`. This can be done by specifying each field and the value you want to change it to as a Name and Value pair, following standard Matlab syntax. This explanation section can also be duplicated within Matlab by publishing or reading the `brainframe_Help.m` script.

- `voxUreg`: Binary flag for per-voxel or per-region visualizations. Human & mouse defaults are both 1. A value of 1 indicates per-region visualizations, and a value of 0 indicates per-voxel visualizations. 

- `brain_atlas`: 3D volumetric atlas with numeric regional IDs per-voxel. Human default is the 86-region Desikan & mouse default is the 213-region AIBS CCF. 

- `bgcolor`: Image background color. Options are 'k', 'w', 'other' (which produces gray). Human & mouse defaults are 'k'.

- `savenclose`: Binary flag indicating whether the user wants images of on axis views saved and the GUI to close (value = 1), or whether the user wants to open the GUI (value = 0). Mouse and human defaults are both 0.

- `img_labels`: Your desired filename, as a string, if savenclose == 1. Default is 'yourfilename'. You can also specify a full filepath.

- `img_format`: Desired image file format as a string, if savenclose == 1. Default is 'png'.

- `data`: Desired data input, specified as either a vector of per-voxel or per-region entries. Per voxel data can also be specified as a 3D matrix. Default for humans is per-region tau PET pathology data and the default for mice is per-region semi-quantitative tau IHC/IF data.

- `nbin`: Number of bins per voxel data is divided into for colormap visualization. This field only applies to per-voxel data. Default is 1.

- `nreg`: Number of regions to go through in specified atlas. This must be a number equal to  the number of unique region IDs in brain_atlas. This field is only relevant for per-region data. Default is 86 for humans and 426 for mice. 

- `region_groups`: This is a vector specifying the group each region is part of for colormap purposes. This vector must be the of length nreg. Default is an 86-element vector of 1s for humans and a 426 element vector of 1s for mice.

- `sphere`: Binary flag specifying whether or not to visualize spheres at region centers. This field is relevant only for per-region visualizations. Default is 0.

- `sphere_npts`: This field specifies how many points are used to construct the spheres, relative to each sphere's radius. Default is 35. 

- `centered`: This is a two element vector, with fields specifying different
elements of pointcloud functionality. The first element is a binary flag
specifying whether or not the pointclouds are centered at each regional
centroid. The second element denotes the degree of centering, with higher
numbers resulting in more centering. Default is [1 2].

- `cmap`: The colormap used in either per-region or per-voxel visualizations. For per-region visualizations, this must be a number of unique region_groups X 3 matrix of RGB vectors per row. For per-voxel visualizations, this must be an nbin X 3 matrix of RGB vectors per row. Default is [1 0 0]. 

- `xfac`: Universal multiplier for sphere radius sizes or point cloud size and density for both per-region and per-voxel visualizations. Default = 1.

- `pointsize`: Specifies the size of points in the visualizations. Default is 50 for humans and 1 for mice. 

- `iscon`: Binary flag specifying whether to visualize connectivity. Default is 0. 

- `conmat`: Connectivity matrix that is nreg X nreg. Default is 86 X 86 Desikan connectome for humans and the 426 X 426 AIBS connectome for mice. 

- `con_rescale`: Universal connectivity multiplier that scales the number of and spread of the ellipses plotted to visually simulate neural connectivity. Ellipses are visualized per region pairs in a number proportional to the C(i,j) entry of the connectome. The default for humans is 1 and the default for mice is 0.01.

- `con_width`: Specifies the line width of each ellipse visualized. The default is 0.01.

- `con_regiongroups`: A vector that specifies region groups as integers for the connectivity visualization colormap. This field works analogously to region_groups. Defaults are the same as for region_groups.

- `con_cmap`: The colormap for connectivity visualizations. This field works analogously to cmap. The default is [0 0.447 0.741].

- `con_arch`: This field specifies the degree of curvature in the ellipses. Default is 0.5.

- `conarrow_WL`: Specifies the width and length of the cone arrows used to indicate direction of connection between each region pair. Defaults are [1.5 2.5] in humans and [0.5 0.8] in mice. 

## 3. Brief Examples
Below are some brief examples on the basics of the functionality of the Brainframe package, including using the `brainframe_inputs_human` and `brainframe_inputs_mouse` functions to create user specified input_struct objects. For more extensive examples through various probable use cases, including more human use cases, please see the `brainframe_Help.m` script. For a breakdown of the default input_struct objects that are created under the hood inside the `brainframe_inputs_human` and `brainframe_inputs_mouse` functions, please refer to the `brainframe_humandefault_creator.m` and `brainframe_mousedefault_creator.m` as example scripts and **default_mouse** and **default_human.mat** as example input_struct files. These default examples can be found within the **ExampleDefaults** folder.

- **Visualizing the human default input_struct**:

	- Change this to alter the path you load from:

	> matpath = cd;

	- Create input_struct:

	> input_struct = brainframe_inputs_human(matpath);
	> brainframe(input_struct)

	![human_default](/ExampleImages/brainframe_Help_01.png)

- **Visualizing the mouse default input_struct**:

	- Change this to alter the path you load from:

	> matpath = cd;

	- Create input_struct:

	> input_struct = brainframe_inputs_mouse(matpath);
	> brainframe(input_struct)
	> view([-1 0 0]);

	![mouse_default](/ExampleImages/brainframe_Help_02.png)

- **Saving and closing images, rather than opening image GUI**:

	- Create input_struct with savenclose set to 1 rather than 0:

	> matpath = cd;
	> input_struct = brainframe_inputs_mouse(matpath,'savenclose',1)
	> brainframe(input_struct)

	- `brainframe.m` will now create 3 on-axis view image files but will not render an image in the Matlab figure GUI. These files will be saved into cd or into the specified filepath.

- **Using the `region_groups` and `cmap` fields to make colorful per-region visualizations**:

	- The below code creates region_groups based on major region IDs:

	> reggroups = zeros(213,1);
	> amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57; 
	> ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; 
	> pal = 150:157; pon = 158:170; str = 171:178; tha = 179:213;
	> reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
	> reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
	> reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
	> reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
	> reggroups(tha) = 13;
	> reggroups = [reggroups;reggroups];

	- Next we create the cmap based on region_groups:

	> cmap = hsv(length(unique(reggroups)));

	- Finally we specify input_struct fields to be customized using a Name and Value pair argument, as is usual in Matlab. This creates a user specified input_struct with the attendant visualization below:

	> matpath = cd;
	> input_struct = brainframe_inputs_mouse(matpath,'region_groups',...
	> reggroups,'cmap',cmap);
	> brainframe(input_struct)
	> view([-1 0 0]);

	![mouse_regioncolors](/ExampleImages/brainframe_Help_05.png)

- **Visualizing only certain regions, rather than all of them**:

	- Note in the below example that regions are visualized as spheres
	- Note in the below example that reggroups and cmap from above are used again here to reset `region_groups` and `cmap`.
	- Note that sphere radii are being reset using `xfac`.
	- Note that this visualizes hippocampal regions only!

	- Set regions you don't want visualized to 0 in the data vector, as below:

	> datavec = zeros(426,1);
	> datavec(27:37) = 1;

	- Reset the fields you want to reset, including using reggroups and cmap from the prior example to reset `region_groups` and `cmap`:

	> input_struct = brainframe_inputs_mouse(matpath,'region_groups',...
	> reggroups,'cmap',cmap,'xfac',0.075,'sphere',1,'data',datavec);
	> brainframe(input_struct)
	> view([-1 0 0]);

	![mouse_selectregions](/ExampleImages/brainframe_Help_11.png)	

- **Visualizing connectivity using the selected regions from above**:

	- Turn on connectivity visualization binary flag:

	> input_struct.iscon = 1;

	- Zeroing out the diagonal:

	> input_struct.conmat = input_struct.conmat - ...
	> diag(diag(input_struct.conmat)); 

	- Zeroing out all non-hippocampal connectivity data:

	> input_struct.conmat(1:26,:) = 0;
	> input_struct.conmat(27:37,[1:26 38:end]) = 0;
	> input_struct.conmat(38:end,:) = 0;

	- This line thresholds conmat:
	
	> input_struct.conmat(input_struct.conmat<0.5*...
	> mean(nonzeros(input_struct.conmat))) = 0;

	- Visualize connectivity:

	> brainframe(input_struct);
	> view([-1 0 0]);

	![mouse_selectregions](/ExampleImages/brainframe_Help_12.png)

- **Per-voxel visualizations**:

	- First load in and select the example per-voxel data:

	> matpath = cd;
	> load([matpath filesep 'PerVox_ExampleData.mat'],'pervoxdata');

	- Here we use Pvalb+ interneuron distributions

	> datinput = pervoxdata.Pvalb;

	- Co-register your data to the reference atlas!

	> datinput = imresize3(datinput,[133 81 115]);
	> datinput(datinput<0) = 0;

	- Create your input_struct and visualize per-voxel data:
	
	> input_struct = brainframe_inputs_mouse(matpath,'voxUreg',
	> 0,'data',datinput,'nbin',5,'cmap',autumn(5));
	> brainframe(input_struct);
	> view([-1 0 0]);

	![mouse_selectregions](/ExampleImages/brainframe_Help_13.png)

- **More examples and further information**: for more examples and further information about various use cases, including some human use case examples, please see the `brainframe_Help.m` file. 

## 4. Best Practices, Citations, & Contact
This section recaps best practices, how to use & cite this package, and how to contact the creator. 

- **Best Practices**: Download the contents of this package and keep all scripts in main folder together. Keep all items in main folder until comfortable with use. Moving items around is possible and can be useful, but is not recommended unless the user is both very familiar with Matlab and has experience with Brainframe. Further help and more use case examples can be found in `brainframe_Help.m` and default examples can be found in the **ExampleDefaults** folder.

- **Citations**: Brainframe can be cited by the Github repository URL and the DOI associated with it. When citing Brainframe, please also cite `arrow3.m`. Citation information and other information about [arrow3.m](https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3) can be found on the Mathworks exchange site in the hyperlinked text.

- **Contact**: to contact the creators of Brainframe, Chris Mezias & Justin Torok, please use the following email addresses:
	-Chris Mezias: cmezias@gmail.com
	-Justin Torok: jlt46@cornell.edu
