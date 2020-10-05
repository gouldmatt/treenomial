## Test environments
* local OS X Catalina install R release
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 2 NOTEs:  
    
* checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      libs   5.8Mb

* checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'treenomial-Ex_i386.Rout' 'treenomial-Ex_x64.Rout'
    'examples_i386' 'examples_x64' 'tests_i386' 'tests_x64'
  
These files/directories were created by R-hub.

## Downstream dependencies
There are currently no downstream dependencies for this package. 
