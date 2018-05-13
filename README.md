volume
======

A command-line tool to work with volume data.

<table>
<tr><td>-i str</td>                                   <td>input volume</td></tr>
<tr><td>-o  str</td>                                  <td>output volume</td></tr>
<tr><td>-selectVolume int</td>                        <td>selects the n-th volume in a nifti file with many volumes	</td></tr>
<tr><td>-largest6con</td>                             <td>largest 6 connected component</td></tr>
<tr><td>-dilate int</td>                              <td>dilate(size)</td></tr>
<tr><td>-erode  int</td>                              <td>erode(size)</td></tr>
<tr><td>-compress   float str</td>                    <td>cosinus transform compress rate coeff_file</td></tr>
<tr><td>-hist   int</td>                              <td>hist(#bins)</td></tr>
<tr><td>-matchHist  str</td>                          <td>matchHistogram(another_mri)</td></tr>
<tr><td>-stats</td>                                   <td>stats, returns mean, std, min, max</td></tr>
<tr><td>-tiff   str str str float</td>                <td>write slice as tiff file. Args: path, cmap, ori {x, y, z}, slice</td></tr>
<tr><td>-info</td>                                    <td>information: dimensions, data type, pixel size</td></tr>
<tr><td>-threshold  float int</td>                    <td>threshold(level,direction)</td></tr>
<tr><td>-volume</td>                                  <td>calculate volume</td></tr>
<tr><td>-new</td>                                     <td>create a new volume with dimx,dimy,dimz,pixx,pixy,pixz,offx,offy,offz</td></tr>
<tr><td>-zigzag</td>                                  <td>print volume values in zigzag order</td></tr>
<tr><td>-decompose  str</td>                          <td>decompose(basename) a volume with many values into volumes with one single value</td></tr>
<tr><td>-strokeMesh str</td>                          <td>set the surface of the mesh (text format) at input path to value=max+1; mesh needs to be in voxel coordinates</td></tr>
<tr><td>-surfaceNets level path</td>                  <td>extract isosurface from the volume at the indicated level using the surface nets algorithm, save at the indicated path</td></tr>
<tr><td> -sampleMesh str1 str2</td>                   <td>sampleMesh(mesh_path, result_path) save the volume values at the vertices of the mesh pointed by the file path to the result path</td></tr>
</table>
