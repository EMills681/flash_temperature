Predict the flash and bulk temperature of a polymer gear (gear ratio: 1:1)

See Python_Abaqus flowchart.png

Default material POM

To use this code in Abaqus, save all scripts to your Abaqus work directory OR in Abaqus set your work directory to the folder where these files are located.

Then, in Abaqus go to file -> run script... and select:

          1) input_flash_temp.py
              to create a model with 10 heat partitions on the meshing surface
              boundary conditions are applied to the heat partitions in one heating step for the whole duration

          2) input_100step.py
              to create a model with 10 heat partitions on the meshing surface
              a new heating step is applied for each partition and 1 boundary condition is applied to each step
              the heating steps run through the 10 boundary conditions along the tooth face and then a rotation step removes the boundary
                conditions whilst the other gear teeth are meshing before the steps ar repeated. 100 steps : 9 cycles
                
                
Enter inputs (use suggested if you want to see how it works). The model only works with sensible inputs
Check model (and mesh) then run a job.
