### Build

To build and test the DEAP FEMA Parallized work you can do the following steps.  

From the DEAP directory run:

```
make
```

This will build the matlab code into a binary that you can run later.

### Generate Cache

Next you will need to build the cache for the workers.  to do this run:

```
./run_gencache.sh
```

This should build the cache files in the directory data_v5/FEMA_cache.  

### Run Workers

Now we can start the workers by running:

```
./start_workers.sh 10
```

The param 10 tells it how many workers to start.  This number should match the "nfrac" param set in the gen_cache.sh and run_client.sh script.  If you want to change the number of workers all three scripts must match.

This will start the workers and block this shell while the workers are running.  You can hit CRTL-c to cause the workers to exit and any time.  

### Run Job

For this step you will need to have a second shell to run the client from.  So in a new shell (with the workers still running in the old shell) run the following:

```
./run_client.sh
```

This should execute the job and test the results.  The outputs will be in the data_v5/output folders.  

### Customizing

Each of the shell scripts has some params that are set at the start.  Care should be taken to ensure that these match.  

The code currently relies on a design matrix, it looks for this design matrix in your home directory as "designmat_v5.txt", you can replace this name in the run_gencache.sh and the run_client.sh scripts.  It is the valriable "fname_design".

To change the phenotype of the analysis you will need to change the start_workers.sh script and the run_client.sh script.  By default these are currently set to RNI/dmri.  This can be changed by changeing the variable "dirname_imaging" and "fstem_imaging" in the run_client.sh scipt and changeing "fstem_imaging" in the start_workers.sh script.  THe run_gencache.sh scripts generates the cache for 3 phenotypes (smri/JA,dmri/RNI,dmri/RND).  When the worker runs it loads the type it was instructed to.  Workers only support one phenotype at a time.  