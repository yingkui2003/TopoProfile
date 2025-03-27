# How to download and use TopoProfile in ArcGIS Pro
The github site includes an ArcGIS Pro toolbox (atbx) file and a python folder, including four python source code files associated with these tools. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://github.com/user-attachments/assets/b6a3676e-53f3-45bb-aa24-c3f9a71e0679)

A zip file of the whole github folder will be downloaded to the local computer. Unzip this file will create a TopoProfiles-main folder with all folders and files. The user can use the toolbox in ArcGIS Pro, check the source codes, and continue improving this toolbox. Note that the source code file has been imported to each tool in the current version, so that the toolbox can be run without the python files. However, if you revise the python files, the tool will not changed until you import the python file to each tool.    

# TopoProfile
TopoProfile is an ArcGIS toolbox to delineate streamlines, generate cross-valley topographic profiles, derive profile integral and a range of other metrics for both one-sided and two-sided profiles. This toolbox includes four tools to: 1) delineate streamlines or valley bottom lines; 2) generate cross sections along the streamlines based on user-specified spacing and width of the cross sections; 3) derive PIs and other related metrics for specified cross sections; and 4) derive the PI and other metrics related to longitudinal profiles. All these tools are written with Python in ArcGIS Pro and can be used for both one-sided and two-sided topographic profiles.

![image](https://github.com/user-attachments/assets/3124be2a-a1a0-49be-8ae1-d06fc4a1bab2)


# Delineate streamlines
The ‘Delineate streamlines’ tool derives a streamline or valley bottom line based on a digital elevation model (DEM) and a specified lowest valley boundary or cross section, which can be defined as an input feature class file (polylines or polygons) or can be digitized on screen. Other required inputs include the minimum source area to start a streamline (km2), the minimum total contribution area of a stream before joining another stream (km2), and the minimum tributary to main valley ratio (Li and Zhao 2022; Li, 2023). The tool also provides three different ways to smooth streamlines: 1) fixed smoothing based on a user-specified smoothing distance; 2) varied smoothing based on Kienholz et al. (2014) for glacier studies; and 3) no smoothing. The output of this tool includes derived streamlines and optional watershed boundaries corresponding to the lowest valley locations. This tool is revised from the PalaeoIce toolbox (Li, 2023).

![image](https://github.com/user-attachments/assets/95a3be4b-d8cb-426d-ae7a-fb4ec2f65da6)


# Generate Cross Sections
The ‘Generate Cross Sections’ tool generates a set of cross sections along streamlines based on a specified spacing and maximum width on each side of the streamline. The inputs of this tool include the DEM, streamline (can be a file or digitized on the screen), the specified spacing along the streamline, the maximum width on each side of the streamline, and the minimum cross-sectional width to save. This tool also includes two optional inputs: one is confining boundaries, such as watershed boundaries, and the other is areas of exclusion such as lakes or glaciers that the user does not want to be crossed by sections. For the streamline network with tributaries, the tool generates the sub-watershed corresponding to each tributary streamline to ensure that cross sections are generated only within the corresponding sub-watershed of the streamline. In addition, this tool provides three options to further constrain cross-section widths: 1) None (no constraint needed); 2) Refine the cross sections by the highest point on each side. This removes sections that potentially cross over the highest point on each side of the valley; and 3) Terminated by a point of convexity on each side. The method to derive convex points along each side of a cross section is based on the algorithm to determine the potential cirque threshold points (convex turning points) along streamline profiles (Li and Zhao, 2022). Please check Li and Zhao (2022) for details. This algorithm produces a list of major convex points on each side of the initial cross section with a fixed width. The cross sections are terminated at the lowest convex points on each side of the valley, above a height threshold. Specifically, a minimum height above the valley floor corresponding to a specified minimum cross section width is required for this method to avoid selection of a convex point close to the valley bottom. The main output is the generated cross-sections along the streamlines. The tool also provides options to divide cross sections in half by streamlines, save the highest/convex points along the cross sections, and export the topographic plot along each cross-section to a specified folder. The end points and topographic plots can be used to manually check, revise, and adjust the cross-sections generated by this method.

![image](https://github.com/user-attachments/assets/819cb69b-a44c-4cfe-83c6-5cd1c2978827)

# Derive Cross Section Metrics
The ‘Derive Cross Section Metrics’ tool derives metrics for cross sections (two-sided profiles). The inputs include the DEM, cross sections (can be a file or digitized on screen), the method to refine cross sections, the minimum cross section width between two sides, and the minimum height above the valley floor when refining cross sections. The outputs include the refined cross sections with derived metrics and optional points to refine the cross sections and half cross-section metrics. Because users may provide their own digitized cross sections, this tool also provides the same three options described in the previous tool to refine the cross sections. To help determine the appropriate cross sections, the tool also generates the highest and convex points on both sides of the valley as an optional output. To help facilitate comparisons of morphologies on each side of the valley, the tool also provides an option to save the half cross-section metrics. Similar to the previous tool, this tool also provides an option to specify a folder in which to save the plots of cross-sectional profiles for further analysis.

![image](https://github.com/user-attachments/assets/162d9473-fa8a-4dee-ad5f-248c3d4b9778)

# Derive One-sided Profile Metrics
The ‘Derive One-sided Profile Metrics’ tool derives the metrics for one-sided profiles, such as river longitudinal profiles, slope profiles, and half valley profiles. The inputs include the DEM, topographic profiles (can be a file or digitized on screen), and whether it is necessary to cut profiles by the highest elevation points. The output consists of the profiles with derived metrics as attributes (Fig. 2d). This tool also provides an option to specify a folder in which to save the plots of the topographic profiles for further analysis.

![image](https://github.com/user-attachments/assets/ecd7a48d-c4e0-4904-8a21-13f0123d6261)


# Cite this work
Li Y., Evans, I.S., Harbor, J., in review. Profile Integral: a robust and flexible metric to measure the concavity/convexity of a topographic profile.

# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu







