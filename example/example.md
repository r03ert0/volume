Example 1. Use the surface nets algorithm to create a mesh from a nii.gz volume

1. Open a .nii.gz volume and display the maximum and minimum values:

````
$ ./volume -i example/six.nii.gz -stats
volume, v4, roberto toro, 9 August 2015
Format: NiftiGZ volume
Tempname: /tmp/nii.g03z4x
mean 46.7929
std 86.1517
min 0
max 253
````

2. Create a mesh at level 10

````
$./volume -i example/six.nii.gz  -surfaceNets 10 ~/Desktop/six.txt
````

3. The mesh is saved in a very simple but non-standard text format. Surface nets creates meshes with a lot of vertices sticking together. To convert the mesh to a standard format and correct the sticking vertices, we use "meshgeometry" (http://github.com/r03ert0/meshgeometry), like this:

````
$ ./meshgeometry_mac -i ~/Desktop/six.txt -fixNonmanifold -o ~/Desktop/six.ply
````

4. Finally, we open this mesh in MeshLab to decimate (Quadratic Edge Collapse Decimation)it and smooth it (Taubin Smoothing). The result is the mesh six-smooth.ply

<script src="https://embed.github.com/view/3d/r03ert0/volume/master/example/six-smooth.stl"></script>
