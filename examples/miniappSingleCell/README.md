# This directory/readme is for use of singleCell example as a miniapp
The reason it exists is for instances where `SENSEI_FOUND` in [CMakeLists.txt](singleCell/CMakeLists.txt) is false (which is the case for SENSEI 4.0.0 at the moment (August 19 2022)) 
If you're having the same issue, this directory provides a way for you to build singleCell as part of the SENSEI build similar to how you'd build the oscillators example. 

1. For starters copy the contents of this directory to `SENSEI/miniapps` - which will replace your original CMakeLists.txt in that directory

2. You need to modify the directory path in [cellFlow.cpp](singleCell/cellFlow.cpp) to match the path to your [cellFlow.xml](singleCell/cellFlow.xml)  [it will be replaced with a variable path down the line but for now since it's a rough draft it is this way]. To do so just search for "tishn" in [cellFlow.cpp](singleCell/cellFlow.cpp) and modify the path in that line to match one your [cellFlow.xml](singleCell/cellFlow.xml) is in. Alternatively you can also see the line in question [here](https://github.com/tish-n/insituBloodFlow/blob/881f0f91558dd26089ec5caf77f22041503153e3/examples/miniappSingleCell/singleCell/cellFlow.cpp#L396) 

3. You DO need to change paths in the [CMakeLists.txt](singleCell/CMakeLists.txt) `SENSEI/miniapps/singleCell/CMakeLists.txt` [they will be replaced with a variable path down the line but for now since it's a rough draft they are this way]

4. Assuming you have a working build of SENSEI 4.0.0 and you've successfully ran the bidirectional oscillator example you should be able to rebuild the SENSEI  

8. Make sure you enable singleCell example by adjusting your cmake options to include `-DENABLE_SINGLECELL=ON`. Alternatively, if your directory heirarchy matches that which is achieved by [setup script](../../setup.sh) feel free to copy the [onlySENSEI.sh](onlySENSEI.sh) script to directory containing your insituBloodFlow directory (once again assuming you used setup.sh to setup everything) and running it from there. 

5. Before you run the executable (which should save to `install/sensei/bin`) check to make sure all of the libraries have been linked properly. To do so run `ldd ./cellFlow` where you have the executable saved. If you have a bunch of libraries that show as `"==> not found"` add them to your library path by running the following commands: 
   `export LD_LIBRARY_PATH=/path/to/missing/library:$LD_LIBRARY_PATH` 

   for me it was: 
   `export LD_LIBRARY_PATH=/home/ntishchenko/scriptTest/insituBloodFlow/install/sensei/lib/python-3.10/site-packages/sensei:$LD_LIBRARY_PATH`
   and
   `export LD_LIBRARY_PATH=/home/ntishchenko/myFork/insituBloodFlow/install/paraview/lib:$LD_LIBRARY_PATH`

6. Once you have that done if you type the following: 
   `echo $LD_LIBRARY_PATH` it should output something like this:
`/home/ntishchenko/scriptTest/insituBloodFlow/install/sensei/lib/python-3.10/site-packages/sensei:/home/ntishchenko/myFork/insituBloodFlow/install/paraview/lib`

7. Keep adding paths to libraries if they're missing. If you're missing one of the core libraries, doublecheck to make sure you have all of the dependencies installed. The terminal output you should get should resemble ones you got when you ran singleCell with palabos/lammps coupling only.

1. Your simulation should be able to initialize and look like the following:  ![image](https://user-images.githubusercontent.com/84354064/185713722-f9aa1a3c-3cbe-47e5-881c-7050e7dc667e.png)


9. Regenerate the `catalystInitializer.py` using ParaView 
  * run the simulation with posthoc enabled to generate a dataset. To do so, change the value in [cellFlow.xml](singleCell/cellFlow.xml) next to posthoc analysis to '1' and next to catalyst analysis to '0'.
  * add a PNG extractor to successfully generate the script 
  * generate a catalyst script in ParaView using the dataset -> load in the dataset -> File -> Save Catalyst State  
  * name it catalystInitializer3.py to be consistent with the [cellFlow.xml](singleCell/cellFlow.xml) file or simply change the name in the [cellFlow.xml](singleCell/cellFlow.xml) file.

Your output should resemble the following ![image](https://user-images.githubusercontent.com/84354064/185712897-aff8744d-5356-4dbc-9ab9-5cad1141ebd1.png)
