This directory/readme is for use of singleCell example as a miniapp
The reason it exists is for instances where SENSEI_FOUND in CMakeLists.txt is false (which is the case for SENSEI 4.0.0 at the moment (August 19 2022)) 
If you're having the same issue, this directory provides a way for you to build singleCell as part of the SENSEI build similar to how you'd build the oscillator example. For now the commit is a rough draft so you might need to change some of the directories in one of the CMakeLists.txt

For starters copy the contents of this directory to SENSEI/miniapps

Now you need to modify the directory path in cellFlow.cpp to match the path to your cellFlow.xml [it will be replaced with a variable path down the line but for now since it's a rough draft it is this way]. To do so just search for "tishn" in cellFlow and modify the path in that line to match one your cellFlow.xml is in

You DO need to change paths in the SENSEI/miniapps/singleCell/CMakeLists.txt [they will be replaced with a variable path down the line but for now since it's a rough draft they are this way]

[Assuming you have a working build of SENSEI 4.0.0 and you've successfully ran the bidirectional oscillator example you should be able to rebuild the SENSEI]  

before you run the executable (which should save to install/sensei/bin) check to make sure all of the libraries have been linked properly. to do so run ldd ./cellFlow where you have the executable saved.

if you have a bunch of libraries that show as "==> not found" add them to your library path by running the following commands: 

export LD_LIBRARY_PATH=/path/to/missing/library:$LD_LIBRARY_PATH

for me it was: 

export LD_LIBRARY_PATH=/home/ntishchenko/scriptTest/insituBloodFlow/install/sensei/lib/python-3.10/site-packages/sensei:$LD_LIBRARY_PATH

and export LD_LIBRARY_PATH=/home/ntishchenko/myFork/insituBloodFlow/install/paraview/lib:$LD_LIBRARY_PATH

once you have that done if you 

echo $LD_LIBRARY_PATH it should output something like this:

/home/ntishchenko/scriptTest/insituBloodFlow/install/sensei/lib/python-3.10/site-packages/sensei:/home/ntishchenko/myFork/insituBloodFlow/install/paraview/lib


