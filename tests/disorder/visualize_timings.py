import matplotlib.pyplot as plt
import numpy as np
#This is the table of times in algortihm old:
#Time for 1 threads:  144071084 microseconds
#Time for 2 threads:  77769634  microseconds
#Time for 3 threads:  54250524  microseconds
#Time for 4 threads:  42241603  microseconds
#Time for 5 threads:  35572385  microseconds
#Time for 6 threads:  32265488  microseconds
#Time for 7 threads:  29077272  microseconds
#Time for 8 threads:  27287518  microseconds
#Time for 9 threads:  25698694  microseconds
#Time for 10 threads: 24265947  microseconds
#Time for 11 threads: 23505353  microseconds
#Time for 12 threads: 23625557  microseconds
#Time for 13 threads: 24480940  microseconds
#Time for 14 threads: 23957628  microseconds
#Time for 15 threads: 23490166  microseconds
#Time for 16 threads: 23147351  microseconds
#Time for 17 threads: 22801197  microseconds
#Time for 18 threads: 22691655  microseconds
#Time for 19 threads: 22345215  microseconds
#Time for 20 threads: 22258808  microseconds
#Time for 21 threads: 22525129  microseconds
#Time for 22 threads: 23288661  microseconds
#Time for 23 threads: 23640912  microseconds
#Time for 24 threads: 25242629  microseconds
#Time for 25 threads: 26527249  microseconds
#Time for 26 threads: 27357247  microseconds
#Time for 27 threads: 27349900  microseconds
#Time for 28 threads: 27513258  microseconds
#Time for 29 threads: 27019550  microseconds
#Time for 30 threads: 27075544  microseconds
#Time for 31 threads: 27204760  microseconds
#Time for 32 threads: 27060650  microseconds
#Time for 33 threads: 26712824  microseconds
#Time for 34 threads: 26599144  microseconds
#Time for 35 threads: 26330746  microseconds
#Time for 36 threads: 26207819  microseconds
#Time for 37 threads: 25932689  microseconds
#Time for 38 threads: 26099687  microseconds
#Time for 39 threads: 25684003  microseconds
#Time for 40 threads: 25919481  microseconds
#Time for 41 threads: 25292227  microseconds
#Time for 42 threads: 25005242  microseconds
#Time for 43 threads: 25106320  microseconds
#Time for 44 threads: 25167409  microseconds
#Time for 45 threads: 25499794  microseconds
#Time for 46 threads: 25456367  microseconds
#Time for 47 threads: 26224751  microseconds
#Time for 48 threads: 26240714  microseconds

#This is the table of times in algortihm new:
#Time for 1 threads:  144894208 microseconds
#Time for 2 threads:  77262801  microseconds
#Time for 3 threads:  53711702  microseconds
#Time for 4 threads:  41966548  microseconds
#Time for 5 threads:  35056128  microseconds
#Time for 6 threads:  30623684  microseconds
#Time for 7 threads:  28841719  microseconds
#Time for 8 threads:  26780783  microseconds
#Time for 9 threads:  25267139  microseconds
#Time for 10 threads: 23678121  microseconds
#Time for 11 threads: 23143896  microseconds
#Time for 12 threads: 23186984  microseconds
#Time for 13 threads: 23453988  microseconds
#Time for 14 threads: 22589325  microseconds
#Time for 15 threads: 22177296  microseconds
#Time for 16 threads: 22507356  microseconds
#Time for 17 threads: 22420556  microseconds
#Time for 18 threads: 22582314  microseconds
#Time for 19 threads: 22318176  microseconds
#Time for 20 threads: 22244671  microseconds
#Time for 21 threads: 22576321  microseconds
#Time for 22 threads: 22861070  microseconds
#Time for 23 threads: 23327785  microseconds
#Time for 24 threads: 24541733  microseconds
#Time for 25 threads: 26163111  microseconds
#Time for 26 threads: 26803355  microseconds
#Time for 27 threads: 26891744  microseconds
#Time for 28 threads: 26678702  microseconds
#Time for 29 threads: 26582280  microseconds
#Time for 30 threads: 26524850  microseconds
#Time for 31 threads: 26637648  microseconds
#Time for 32 threads: 26425435  microseconds
#Time for 33 threads: 26368898  microseconds
#Time for 34 threads: 26339228  microseconds
#Time for 35 threads: 26049208  microseconds
#Time for 36 threads: 25773944  microseconds
#Time for 37 threads: 25564022  microseconds
#Time for 38 threads: 25141991  microseconds
#Time for 39 threads: 24891307  microseconds
#Time for 40 threads: 24667356  microseconds
#Time for 41 threads: 24332523  microseconds
#Time for 42 threads: 24129395  microseconds
#Time for 43 threads: 24021825  microseconds
#Time for 44 threads: 23941234  microseconds
#Time for 45 threads: 23737689  microseconds
#Time for 46 threads: 23527481  microseconds
#Time for 47 threads: 23690486  microseconds
#Time for 48 threads: 24136076  microseconds

#plot the following 2 tables into a figure as a scatter plot where we plot time in microseconds against the number of threads used:
nthreads = np.array(range(1,49))
times_old = np.array([144071084, 77769634, 54250524, 42241603, 35572385, 32265488, 29077272, 27287518, 25698694, 24265947, 23505353, 23625557, 24480940, 23957628, 23490166, 23147351, 22801197, 22691655, 22345215, 22258808, 22525129, 23288661, 23640912, 25242629, 26527249, 27357247, 27349900, 27513258, 27019550, 27075544, 27204760, 27060650, 26712824, 26599144, 26330746, 26207819, 25932689, 26099687, 25684003, 25919481, 25292227, 25005242, 25106320, 25167409, 25499794, 25456367, 26224751, 26240714])
times_new = np.array([144894208,77262801 ,53711702 ,41966548 ,35056128 ,30623684 ,28841719 ,26780783 ,25267139 ,23678121 ,23143896 ,23186984 ,23453988 ,22589325 ,22177296 ,22507356 ,22420556 ,22582314 ,22318176 ,22244671 ,22576321 ,22861070 ,23327785 ,24541733 ,26163111 ,26803355 ,26891744 ,26678702 ,26582280 ,26524850 ,26637648 ,26425435 ,26368898 ,26339228 ,26049208 ,25773944 ,25564022 ,25141991 ,24891307 ,24667356 ,24332523 ,24129395 ,24021825 ,23941234 ,23737689 ,23527481, 23690486, 24136076])

plt.plot(nthreads, times_old, label='old')
plt.plot(nthreads, times_new, label='new')
plt.xlabel('Number of threads')
plt.ylabel('Time in microseconds')
plt.legend()
plt.show()

#calcualte the speedup from these timings and plot:
speedup_old = np.zeros(48)
speedup_new = np.zeros(48)
for n in range(1,49):
    speedup_old[n-1] = times_old[0]/times_old[n-1]
    speedup_new[n-1] = times_new[0]/times_new[n-1]
    
plt.plot(nthreads, speedup_old, label='old')
plt.plot(nthreads, speedup_new, label='new')
plt.xlabel('Number of threads')
plt.ylabel('Speedup')
plt.legend()
plt.show()