# @ job_type = parallel
# @ job_name = testbugs
# @ node = 4
# @ tasks_per_node = 16
# @ wall_clock_limit = 1:00:00
# @ notification = never
# @ input = /dev/null
# @ output = out.$(jobid)
# @ error  = err.$(jobid)
# @ queue

export OMP_NUM_THREADS=2
cd $LOADL_STEP_INITDIR
./test.par
