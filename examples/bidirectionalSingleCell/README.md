# HOW TO RUN THIS EXAMPLE

To run this example, do the following steps:

if you haven't already, clone this repo using git clone 
`git clone https://github.com/tish-n/insituBloodFlow`
then 
`cd insituBloodFlow`
then
`git checkout devtest`

now modify the cellFlow.cpp file with appropriate path for xml config file. It should point to the singleCell4.0/cellFlow.xml

modify the CMakeLists.txt in singleCell4.0 directory with the path to your LAMMPS installation. Ideally, if you use the setup.sh script it will point to the directory you cloned this repo to
`path/to/where/you/cloned/this/repo/to/insituBloodFlow/src/lammps/src`

now just run the [setup.sh](../../setup.sh) script from the directory you cloned
