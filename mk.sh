#! /bin/sh
make
echo "\nWorking with $1\n"
echo "======= Sequential ======="
./bsp_mcl_seq $1
echo 
for p in 2 4
do
	echo "======= Parallel p=$p ======="
	./bsp_mcl $p $1
	echo 
done