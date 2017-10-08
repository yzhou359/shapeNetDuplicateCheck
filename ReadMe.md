
# ShapeNetDuplicateCheck Project Overview

## C++ implementation of checking duplicates in ShapeNetCore v2 dataset.

1. Create poison-sample bash file & re-computer vertex norm (optional) bash file.

2. Encode 3D point clouds to 6x6x6 grid representation.

3. Casually find K (=100) most similar neighbors in the same category.

4. Accuratly find distance between neighbors (point-to-plane distance w.r.t. bounding box size) and judge duplicates (dis<0.01%).

5. Find total number of duplicates based on duplicate-xxx file.


## Out of this solution:

1. create list.txt in ShapeNetCore.v2 each folder [before step 1]

	\> ~\ShapeNetCore.v2\ShapeNetCore.v2\$categoryNum\ dir /b | findstr ply > list.txt

2. create render images for ShapeNet Models [before step 1]

	\> ~\ShapeNetCore.v2\lun_script\render run-render.bat categoryNum

3. create distance folder for each category [before step 4]

	\> ~\ShapeNetCore.v2\ShapeNetCore.sample\$categoryNum\ mkdir dis

4. copy list.txt & hash.txt file [before step 4]

	\> ~\ShapeNetCore.v2\ShapeNetCore.sample\$categoryNum\ cp 'list.txt' 'list - Copy.txt' | cp 'hash.txt' 'hash - Copy.txt'

5. create list.txt in dis folder [before step 5]

	\> ~\ShapeNetCore.v2\ShapeNetCore.sample\$categoryNum\dis\ dir /b | findstr duplicate- > list.txt



## Visualization Website: 

Check it [here](http://antares.cs.umass.edu/project_data/AdversarialMonsters/Duplicate/ShapeNetDuplicate/index-bsr-exterior.php?category=02691156).

