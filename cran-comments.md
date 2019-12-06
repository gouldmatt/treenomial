## Test environments
* R-hub solaris-x86-patched (r-patched)
* R-hub ubuntu-gcc-release (r-release)
* R-hub fedora-clang-devel (r-devel)
* R-hub windows-x86_64-devel (r-devel)
* R-hub fedora-gcc-devel (r-devel)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 3 NOTEs:  

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Matthew Gould <mgould@sfu.ca>'

  Days since last update: 0
  
Changed all input checks to use is() rather than checking with class() fixing errors in the R 4.0.0 release where matrix objects now also inherit from class "array". This update was done to retain the package on CRAN after being notified of the problem with the treenomial 1.0.1 version. 
  
  Possibly mis-spelled words in DESCRIPTION:
    Liu (19:35)
    phylo (20:29)
    
These words are spelled correctly. 
    
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
