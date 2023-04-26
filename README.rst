==================
MicroMPN: A Python command line program for automating most probable number (MPN) estimates in microtiter plates 
==================

Authors
------------------
Adam Rivers, United States Department of Agriculture, Agrcultural Research Service
Karla Franco Meléndez, United States Department of Agriculture, Agrcultural Research Service


About MPN
------------------
Most probable number (MPN) is a  way of estimating now many viable bacteria are in a sample. 
A series of dilutions are made with replication and the  number of  tubes  or wells with bacterial growth is recored.
This presence/absence and dilution information is used to calcualte the mos probable number of bacteria in a starting sample.

With liquid handling robots and plate reader thsis method can be baster than conventniona plate counting or 
spotting for some applications. 

The addition of a florescent mmarker gene into a microbe of interest can be 
used to screen for conditions that alter the growth of a specifi microbe in a complx community 
even when selective media does not work.

About MicroMPN
---------------


This python package has a command line program processing data from a plate reader. Users provide two types of files: 
A plate layout file specified in Wellmap format and a  CSV plate reader file with columns for the Plate, 
the well or row and column and the optical value used for determining growth.


Developers of Python Packages 
who need an MPN calculation may want to use MicroMPN's mpn estimation function.
This is to my knowledge the only python implmentation of MPN estimation in the Python Package Index.
This function is derived from Martine Ferguson and John Ihrie's MPN R package:

Martine Ferguson and John Ihrie. (2019). MPN: Most Probable Number and 
Other Microbial Enumeration Techniques. R package version 0.3.0. https://CRAN.R-project.org/package=MPN .

Installing
---------------

The package can be installed from the python package index 

..  code-block:: bash
    pip install micrompn


Getting Started
---------------

