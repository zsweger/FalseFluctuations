`setup 64bits`
`make` #compile
`make clone` # generate jobs 
`nohup ./queue_daemon > nohup.out &`  # submit jobs 

`conf.mk`: variable `NUM` means number of jobs
`urqmd.cc`: line 46 - 75 are configuration for urqmd
`queue_daemon`: qnum means maxmium jobs running number
