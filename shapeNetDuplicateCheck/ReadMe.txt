========================================================================
    CONSOLE APPLICATION : shapeNetDuplicateCheck Project Overview
========================================================================

C++ implementation of checking duplicates in ShapeNetCore v2 dataset.

1. Create poison-sample bash file & re-computer vertex norm bash file.


Out of this solution:

1. create list.txt in ShapeNetCore.v2 each folder [before step 1]
	>> ~\ShapeNetCore.v2\ShapeNetCore.v2\$categoryNum\ dir /b | findstr ply > list.txt

2. create render images for ShapeNet Models [before step 1]
	>> ~\ShapeNetCore.v2\lun_script\render run-render.bat categoryNum

3. create distance folder for each category [before step 4]
	>> ~\ShapeNetCore.v2\ShapeNetCore.sample\$categoryNum\ mkdir dis

4. copy list.txt & hash.txt file [before step 4]
	>> ~\ShapeNetCore.v2\ShapeNetCore.sample\$categoryNum\ cp 'list.txt' 'list - Copy.txt' | cp 'hash.txt' 'hash - Copy.txt'




Visualization Website:
	http://antares.cs.umass.edu/project_data/AdversarialMonsters/UserStudy/ShapeNetDuplicate/index.php
/////////////////////////////////////////////////////////////////////////////
