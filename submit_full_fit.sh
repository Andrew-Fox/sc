#!/bin/sh
# -N paramFit$1

# The job should be placed into the queue 'all.q'.
# -q all.q

#$ -cwd

# Redirect output stream to this file.
# -o MCMCSolutions/paramFit$1.dat

# Redirect error stream to this file.

#$ -e mcmc_error.dat

# The batchsystem should use the current directory as working directory.
# Both files will be placed in the current
# directory. The batchsystem assumes to find the executable in this directory.
# -cwd

# request Bourne shell as shell for job.
#$ -S /bin/sh

# print date and time
date

# spython is the server's version of Python 2.5. Using python instead of spython causes the program to run in python 2.3
/usr/local/bin/matlab -nodisplay -nosplash -r "test $1" 
# print date and time again
date
