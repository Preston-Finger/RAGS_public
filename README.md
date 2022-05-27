# RAGS
## Rapid Assessment of Galaxy Spectra

## Setup:

1. Click green code button. 

2. Download zip. 

3. Move the downloaded zip file to where you want to work 

4. unzip file by doing the command unzip "folder name"

## Before Classifying

1. Grab these files below from these links:

Pears South Spectra Tarball:
 
https://mnscu-my.sharepoint.com/:u:/g/personal/vd3712rv_minnstate_edu/EcDHQoD1jRhKgwO3vPVOh5QBCAWRPcM823w68LcDkRlCKA?e=7opo2x

Goods South mosaic:
 
https://mnscu-my.sharepoint.com/:u:/g/personal/vd3712rv_minnstate_edu/ERrASBQm1ORLhWI3tIst-zgBCuM9f9OagtaEngtxuv0ezA?e=PsDm0W


2. Move Pears South Spectra Tarball from downloads folder to h_pears_total_south in the RAGS_public_main folder

             and 

   Move hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits to master_files folder in RAGS_public_main folder


3. On mac open a terminal and navigate to the h_pears_total_south folder and type ls to confirm pears_south_tar.gz is there

4. enter in the command tar -xvzf h_pears_total_south_tar.gz

5. The file will unzipped and the .fits files will be there. 

6. Repeat these steps with the pears_north_tar.gz

7. Move hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits to master_files folder in RAGS_public_main folder



## Classifying

1. From terminal in the RAGS_public_main folder type the command python RAGS.py  NOTE: or however you run python on your mac

2. Enter in your initials or some sort of designator. 

3. Click or hit enter to begin inspecting objects. 

4. Follow the directions provided for how to classify galaxy spectra.

5. If needed to quit and continue later hit escape and then quit on the window. Upon relaunching program it will continue where last left off of.

6. When program is complete with the list of objects to classify a done message will appear. Hit enter and the program will close.


## After Classifying

1. In the users folder you will see south_Classification_XXX.txt (or north depending on which set classifying) and the XXX being the initials of the user. 

2. This file has the ObjectID, RA, DEC, z_match, Imag, Imag_error, and Classification you assigned. 

3. Send this file to --> ________________


Thank you
