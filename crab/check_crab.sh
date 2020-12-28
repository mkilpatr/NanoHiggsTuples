#!/bin/bash

mypath=crab_projects_mc
for i in ${mypath}/*; do
	echo $i
	crab status $i
done
