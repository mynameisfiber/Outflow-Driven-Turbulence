#!/bin/bash

counter=1
while [ /dev/true ]
do
	mkdir "run-`printf "%.8d" ${counter}`"
	cd "run-`printf "%.8d" ${counter}`"
	cp ../isoth3d .
	./isoth3d | cat > log
	echo "====================" >> ../log
	echo "Run ${counter}" >> ../log
	echo "--------------------" >> ../log
	tail -f 15 log >> ../log
	rm -rf isoth3d
	cd ..
done