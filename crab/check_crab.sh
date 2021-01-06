#!/bin/bash

mypath=crab_projects_mc_v2
for i in ${mypath}/*; do
	echo $i
	crab status $i
done
