import os
#TEGAN: eventually put here a proper function to run the tests based on input arguments
#for now, just comment out the ones I don't want to do

#ideally I should also do a check to see if any of the source files or input files have change
#since the image was last produced and only run  if the have changed. 

os.system("./run_iso.sh")
#os.system("./run_gr_iso.sh")
#os.system("./run_thomson_iso.sh")
#os.system("./run_compton_iso.sh")