## Test environments
* local OS X Mojave install, R 3.6.1
* R-hub macos-elcapitan-release (r-release)
* R-hub windows-x86_64-devel (r-devel)
* R-hub ubuntu-gcc-release (r-release)
* R-hub fedora-clang-devel (r-devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 3 NOTEs:

❯ On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Matthew Gould <mgould@sfu.ca>'
  
  New submission

❯ On windows-x86_64-devel (r-devel)
  checking for non-standard things in the check directory ... NOTE
  Found the following files/directories:
    'treenomial-Ex_i386.Rout' 'treenomial-Ex_x64.Rout'
    'examples_i386' 'examples_x64' 'tests_i386' 'tests_x64'

These files/directories were created by R-hub. 

❯ On ubuntu-gcc-release (r-release)
  checking installed package size ... NOTE
    installed size is  5.9Mb
    sub-directories of 1Mb or more:
      libs   5.8Mb

## Downstream dependencies
There are currently no downstream dependencies for this package. 
