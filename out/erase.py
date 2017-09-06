#!/bin/sh

a = 'rm submit_full_fit_'

import subprocess
import time



for i in range(0,100):
	b = str(i)
	c = a+b+'.error'
	print c
	# Send the command in a shell    
        subprocess.call( c, shell=True )
        # Wait a couple seconds
        time.sleep(.2)





