To build and test the DEAP FEMA Parallized work you can do the following steps.  

From the DEAP directory run:

```
make
```

This will build the matlab code into a binary that you can run later.

Next you will need to build the cache for the workers.  to do this run:

```
./run_gencache.sh
```

This should build the cache files in the directory XXX.  

Now we can start the workers by running:

```
./start_workers.sh 10
```

The param 10 tells it how many workers to start.  This number should match the "nfrac" param set in the gen_cache.sh and run_client.sh script.  If you want to change the number of workers all three scripts must match.

This will start the workers and block this shell while the workers are running.  You can hit CRTL-c to cause the workers to exit and any time.  This is where you need to have a second shell to run the client from.  So in a new shell (with the workers still running in the old shell) run the following:

```
./run_client.sh
```

This should execute the job and test the results.  The outputs will be in the XXX folders.  

Each of the shell scripts has some params that are set at the start.  Care should be taken to ensure that these match.  