
Matclust Help File

	Matclust is a GUI program that runs in Matlab and allows the user to draw polygons or boxes around multidimensional data points.  This 'clustering' of data is used often for multi-electrode recordings of neural data. The program has been tested on Windows, Linux, and Mac platforms.

Features:

	-Manual 2D clustering with polygons
	-Seeded semi-automated clustering
	-3D rotation of feature space
	-Filtering of points to precisely define each polygon?s scope
	-Easy to use interface (context menus, undo states, zooming and panning)

To Install

Put the Matclust folder anywhere you like on your local harddrive.  In matlab, change the directory to the Matclust folder and type: 

matclustsetup

This will compile all the needed files.  If you get an error saying that you don't have a compiler set, you'll need to go to the matlab website and follow their directions to download a free C compiler. Then, add the Matclust folder to you matlab path.  
	

Starting Matclust

	Once Matclust is installed on the computer, simply type 'matclust' in the Matlab prompt to start the program.  You can also feed Matclust an input variable by typing 'matclust(yourvariable)'.  The input variable can be of a few different formats:


1) a matrix.  If the input is simply a matrix, then matclust will consider each row a data point, and each column represents the different parameter to cluster.  Matclust will name each parameter Column 1, Column2, ...
2) a structure.  Here, one of the fields must be called 'params', which contains the matrix described above.  You can also include a few other fields: 

paramnames  -- an N by 1 cell array, with N equal to the number columns in the data matrix.  Each cell contains a string with the name of the corresponding parameter.

filename ? this is a string containing the name of a file containing any original data. User written helper programs can use this info if needed.

customvar ? this field can contain any data that user-written helper programs can use. This data will be stored in the global variable 'clustdata' as clustdata.customvar.


Opening files from Matclust

	To open a raw data file saved to disk, go to File->Open->Raw parameter file. This will allow you to browse to the wanted filename.  The file should have one variable stored inside with the name 'filedata'.  This variable can take on any of the forms described above.
	Once data has been saved from matclust as a matclust file, it can be opened again using File->Open->MatClust file.
	The user can open up multiple files at the same time.  Switching between the files is accomplished using the 'File->Currently open files' menu.  Switching between files may be slow because the backgrounded files are temporarily written to disk.

Drawing clusters

	To create a cluster, you must first select which cluster number to use.  Press one of the colored buttons on the right control panel to select the active cluster.	
      Clusters are created by using either the square drawing tool, the polygon drawing tool, or the ?magic wand? tool (all located above the graph window).  With the square tool, simply click, drag and release to create a box around the desired points.  With the polygon tool click at each vertex, and double click the final point to release the polygon.  You can also use shift+click or the enter key to release the polygon.  If you want to cancel a polygon while in draw mode, press the escape key. With the magic wand, click on the center of a cluster and a polygon will be drawn for you.  This tool works best for clusters that are well-separated from other clusters. For clusters that are not well-separated, drawing the polygons manually is recommended. To the right of the magic wand button, a slider control is used to tell the wand program how many points are in the cluster (from few to many).
	When a new cluster is created, the area next to its colored button changes from grey to the cluster's color.  Two radio buttons also appear inside this area.  These buttons allow the user to quickly manipulate how the cluster is visualized on the screen.  Normally, when both buttons are not activated, the cluster points are colored the same color as the cluster box, and the cluster boxes are visible.  When the left radio button (or the 'Hide' button) is activated, the points for that cluster disappear.  This can make it easier for the user to see other clusters that occupy the same space.  The right button (the 'Exclude' button) makes it seem like the cluster was never drawn by turning the points back to the cluster 0 color, and making the cluster boxes invisible.  When the cluster 0 is made hidden (by toggling the Cluster 0 button) the excluded cluster disappear along with the rest of cluster 0.  To the right of the cluster controls, the user also has the option to hide or exclude all of the clusters.  This action can be reversed by pressing the 'Release' button.
	Underneath the 'hide all' and 'exclude all' buttons, a list box displays where all of the cluster boxes are located for the currently active cluster.  When a new box or polygon is drawn, a new item will appear in this list.  This item will look something like: 'x(2), y(3), t(1), o(1 3 4)'.  This example item says that the box was drawn when x was set to parameter 2, y set to parameter 3, time filter set to number 1, and that filters 1, 3, and 4 were activated in the 'other filters'.  Clicking on this item will take the user back to that configuration to view the box.
	The user can also manipulate drawn clusters using the arrow tool (located above the display window).  Use this tool to click on any cluster box.  The box or polygon will then highlight.  The polygon can then be dragged to any location.  The user can also click and drag individual vertices of a polygon.  If the DELETE key is pressed when a polygon is highlighted, the polygon is deleted (on a Mac, it?s command+DELETE).  This action will not delete the other polygons belonging to the cluster.      

Cluster Options

	The user can right-click on the colored cluster ID button for a particular cluster and a context menu will appear with a selection of commands that can be performed on the cluster.  Here, you can copy the cluster to another cluster number, delete all or parts of the cluster, change the cluster's color and display layer order, or perform any custom-designed operations located in the 'Tools' submenu.  


Making and using time filters

	To make use of time filters, the first column of the data matrix must be the time of each data point.  To create a new time filter, click on an empty number in the ?Time? list.  When the number becomes highlighted, right-click on it and a menu will appear.  Choose 'Add time filter' and a window will appear allowing the user to enter the name, start time, and end time for the filter.  The times must be entered in hours:minutes:seconds format.  Click 'ok' and a new filter should appear in the time filters list.  If you double-click on this filter, the small double arrows in front of 'All points' will move to the new filter, and only the points within the entered times will be shown.  Any clusters created when the filter is on will exclude all points that are excluded in the filter.  When the radio button above the time list is on, filters in the 'Time' list can only be activated one at a time.  However, the user can de-activate this button (also in the ?Options? menu) and multiple times can be viewed simultaneously (combined with OR logic).  In this mode the user can double-click on and off each time filter independently.  The user can also draw polygons when multiple times are selected?this will copy the drawn polygon into each selected time.  Warning: using the multi-time mode can be a bit confusing because multiple polygons from the same cluster can be visible at the same time.  

Making and using other filters

	The 'Filters' list functions similarly to the ?Time? list, except that multiple filters can be combined together, and that the filters are defined though user-written matlab programs.  Unlike in the ?Time? list, when multiple filters are selected, they are combined with AND logic (time filters are combined with OR logic). These filters can be used to cluster subsets of the data defined by any rule the user wants.  Creating a filter in matclust is done by running a matlab program located in the Filters folder in the user's matclust directory.  Running a filter program from matclust is accomplished by clicking on an empty number in the 'other filters' list, then right clicking to make the menu appear.  Clicking 'Add filter' brings up a window that allows the user to name the new filter and browse for the location of the program.
	Writing a filter program requires some knowledge of how data is organized in matclust.  All of matclust's data is organized in four global variables:  clustdata, graphattrib, clustattrib, and figattrib.  Most important of these is clustdata, which stores the data matrix in clustdata.params.  This filter program do not take inputs, instead they access data through these global variables.  The output of the filter should be an N by 1 logical array, where N is the number of rows in clustdata.params.  Here is a simple example of a filter program:

function output = myfilter

global clustdata
output = (clustdata.params(:,3) >= 50);

This will filter out any points that are less than 50 in the third column (points corresponding to 0's in the filter output get filtered out).
	Activating a filter works the same as in the 'Time' box-- simply double click on the filter and double arrows will appear before the name.  Double-clicking again will de-activate the filter.

Using 3D Rotation

	Because manual clustering requires multidimensional data to be projected onto a 2-dimensional space, it can often be advantageous to rotate the view in 3 dimensions to optimize the projection where the polygons are drawn.  To do this, the user can add rotations to the parameter selection box.  Simply right-click on the parameter selection box and choose ?Add rotation?.  A new window will appear that will allow the user to choose three parameters to view in 3D space.  All filters currently active in the main window will also be active in the rotation window.  Use the sliders at the bottom of the window to control the viewing angle.  When a desired angle is found, click the ?Add? button to add the projection to the parameter selection box.

Designing tools

	Just like the filters, the user can write programs that interface with matclust.  Because the matclust variables are stored as global variables, it is possible to change anything about a currently running matclust window.  For general tools that are not specific to a particular cluster number, save the programs in the Tools folder in the matclust directory.  These programs will appear in the 'Tools' drop-down menu at the next startup of matclust.  For programs that perform operations to specific clusters, save them in the ClustTools folder.  These programs will appear in the 'Tools' menu inside the cluster context menu.  The only difference between these two types of tools is that general tool programs do not take any inputs, while the cluster tool programs take one scalar input-- the cluster number.  Both programs should not return any outputs (they have to interface through the global variables). 

Controls

Left,right:		flips through the tool selection buttons
Control+Left, right:	flips through the parameters to plot
Up, down:		flips though the time filters
Shift + hold:		when a polygon is highlighted, it is copied. 
Release shift:		paste the polygon
Escape key:		abort drawing polygon
Delete key:		delete highlighted polygon
0-10 keys:		toggle the ?display? mode of the corresponding cluster
