# Interaction avoidance analysis

The analysis framework presented here is just for review purposes.

## Requirements

- [Matlab](https://www.mathworks.com/products/matlab.html)
- [Matlab Computer Vision Toolbox](https://www.mathworks.com/products/computer-vision.html)
- [Matlab Image Processing Toolbox](https://www.mathworks.com/products/image.html)
- [Python](https://www.python.org/)

## Work flow

1. compute the neighborhood matrix for each cell: `run_Save_Neighborhood_Matrices.m`
2. run permutation: `histoCAT/perm/permutation.py` 
3. pull interaction/avoidance scores: `run_Neighborhood_Master_ours.m`
4. run t-test if multiple groups: `run_Neighborhood_Master_ours_Ttest.m`

## Config Files

To be able to run the Matlab scripts, the user should specify the following parameters in the config file (within the `ConfigFiles` folder. The following parameters need to be included in the file.

1. `input_pano_folder`: the folder where the raw images are located at.
2. `input_cell_segmentation_folder`: the folder where the cell segmentation results are located at. This folder will have the same file organization (structure) as the input, given that the cell segmentation is done using our pipeline. For the cell segmentation pipeline, please refer to the GitHub repository located here: [https://github.com/cellfactory411/Cell-Segmentation](https://github.com/cellfactory411/Cell-Segmentation)
3. `input_celltype_folder`: the folder where the cell pheon-typing results are located at. This folder will have the same file organization (structure) as the input, given that the cell segmentation is done using our pipeline. For the cell segmentation pipeline, please refer to the GitHub repository located here: [https://github.com/walsh-quail-labs/Cell-Phenotyping-Lung](https://github.com/walsh-quail-labs/Cell-Phenotyping-Lung)
4. `input_keyPatientID`: this is the address for a spreadsheet that has keys and patient IDs. The spreadsheet should contain two columns `FileName` and `UniqueImageID`. Generally, we have a one-to-one mapping between filenames and unique image ids but there could exist cases where we have multiple cores of one unique image id and this is predicted by our software and handled at the deeper level of the software.
5. `group_PatientID`: this is the address for a spreadsheet that the groups and patient ids relationships. one other input to our code is that we can image multiple groups for our analysis. These groups can be coming from biological experimental designs and be easily input into the code. If there are no specific groups for your analysis, one can create just a single group with all the patient IDs in the spreadsheet.