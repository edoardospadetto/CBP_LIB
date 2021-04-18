#!/bin/bash

#list test folder and launch executable

for test in $(ls tests | grep .bin)
do 	
	type=${test%??????????.???}
	type=${type^^}
	mpirun -np 4 ./bin/test${type:30} ${test} 5 4

done
