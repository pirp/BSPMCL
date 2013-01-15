all:
	gcc -o bsp_mcl bsp_mcl.c -IMondriaan3.11/src/include -LMondriaan3.11/src/lib -lMondriaan3 -lm