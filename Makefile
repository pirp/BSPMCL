all: par seq
par:
	gcc -o bsp_mcl bsp_mcl.c -IMondriaan3.11/src/include -LMondriaan3.11/src/lib -lMondriaan3 -lm
seq:
	gcc -o bsp_mcl_seq bsp_mcl_seq.c -IMondriaan3.11/src/include -LMondriaan3.11/src/lib -lMondriaan3 -lm