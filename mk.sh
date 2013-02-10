#! /bin/sh
make par
if [ -n "$1" ]
then
	./bsp_mcl $1 $2
else
	./bsp_mcl 2 yeast
fi
