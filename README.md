# NoBA
R-programs for the Nordic and Barents Seas Atlantis
NB!!! the forcing files requires these two steps after being created:
all files (temp, saln, hydro, ice):
ncrename -d t1,t inputfile.nc outputfile.nc

hydrofiles only:
to change dest_k from 1 - 7 to 0 - 6 you'll need to run the following command:
ncap2  -O -s 'dest_k-=1' inputfile.nc  outputfile.nc
