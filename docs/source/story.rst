##############
Why PyEarth?
##############

More than 10 years ago, I relied on both ArcGIS Arc Engine/Object and MATLAB/IDL for various GIS tasks. 

The ESRI ArcGIS system is a powerful GIS platform, but it is not free to use. 
More importantly, it is not easy to run ArcGIS tools in a Linux environment to take the advantage of the parallel computing power of a Linux cluster.

Although MATLAB may be one step ahead of ArcGIS, it is also not free to use. And I had issues with its scaling capability. 
What's more, MATLAB is not designed to work with GIS data. For example, it is not an easy task to mosaic multiple remote sensing images in MATLAB.

IDL was a clear winner at that time. It supports GIS/RS, it can be run in a Linux environment, and it is easy to scale.
By the time, the popular IDL Coyote library (http://www.idlcoyote.com/) is my go-to library for GIS/RS tasks. I even built a library on top of the Coyote library.

IDL solves most of my problems, but it has some limitations other than the license fee. 
Namely, IDL is mostly conceived as functional programming language, which means it is often not used to build large-scale applications. Meanwhile, I started using Python for other tasks. Linking Python with IDL is possiable but trivial. Since I use HPC, that means I need to use two commands to finish a task instead of one, same for debugging.

Besides, I started to work in the FAIR (fairness, accountability, and transparency) domain once I started using Python. 

This is the point I decided to rewrite my IDL library in Python and the birth of PyEarth.
PyEarth is still functional-oriented, meaning, each function is designed to do one thing. PyEarth functions are implemented to support real-world tasks such as reading a raster file.

PyEarth becomes the suporting library for my other projects, such as the PyFlowline project. The principle is very straightforward: if a feature is needed more than once, it should be abstracted as a function in PyEarth.

PyEarth uses the Python ecosystem to do lots of simple things easier. For example, I don't want to cope and paste a large chuck of code to plot a time series data. In PyEarth, that is one single function call with various customization options.

Thank you.