#! /bin/bash
make
echo -e "\nWorking with $1\n"
echo "======= Sequential ======="
for i in {0..9}
do
	./bsp_mcl_seq $1
done
echo 
for p in 2 4
do
	echo "======= Parallel p=$p ======="
	for i in {0..9}
	do
		./bsp_mcl $p $1
	done
	echo 
done