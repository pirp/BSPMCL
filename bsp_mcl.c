#include "BSPedupack1.01/bspedupack.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Mondriaan.h>

int P;

void print_matrix(int s,struct sparsematrix matrix){
	int k;
	for(k=0;k<matrix.NrNzElts;k++){
		printf("%d: (%ld,%ld)=%f\n",s,matrix.i[k],matrix.j[k],matrix.ReValue[k]);
	}
}

void print_own_part(int s, struct sparsematrix matrix){
	int k;
	int start = matrix.Pstart[s];
	int end = matrix.Pstart[s+1];
	for(k=start;k<end;k++) printf("%d: (%ld,%ld)=%f\n",s,matrix.i[k],matrix.j[k],matrix.ReValue[k]);
}

void bsp_mcl(){
	int p,s,t;
	int i;

	bsp_begin(P);

	p= bsp_nprocs();
    s= bsp_pid();
    
    FILE *File;
    struct sparsematrix matrix;
	if (!(File = fopen("test_matrix.mtx", "r"))) printf("Unable to open arc130!\n");
	if (!MMReadSparseMatrix(File, &matrix)) printf("Unable to read arc130!\n");
	fclose(File);

	long length = matrix.NrNzElts;


	long *MatrixI;
	MatrixI = vecallocl(length);

	for(i=0;i<length;i++) MatrixI[i] = 40;

	if(s==1)for(i=0;i<length;i++) printf("%d: %d=%ld\n",s,i,matrix.i[i]);

	bsp_push_reg(MatrixI,length*sizeof(long));
	//bsp_push_reg(matrix.j,length*sizeof(long));
	//bsp_push_reg(matrix.ReValue,length*sizeof(long));
	//bsp_push_reg(matrix.Pstart,(p+1)*sizeof(long));

	bsp_sync();

	if(s==0){
		struct opts Options;

		if (!SetOptionsFromFile(&Options, "Mondriaan.defaults")) printf("Unable to set options from disk!\n");
		if (!ApplyOptions(&Options)) printf("Invalid options!\n");

		/* Distribute the matrix over two processors with an allowed imbalance of 3% and the options provided above. */
		if (!DistributeMatrixMondriaan(&matrix, 2, 0.03, &Options, NULL)) printf("Unable to distribute arc130!\n");
		matrix.i[length-1] = 50;

		for(t=0;t<p;t++){
		bsp_put(t,matrix.i,MatrixI,0,length*SZLONG);			
		//bsp_put(t,&matrix.j,matrix.j,length*sizeof(long),0);			
		//bsp_put(t,&matrix.ReValue,matrix.ReValue,length*sizeof(long),0);			
		//bsp_put(t,&matrix.Pstart,matrix.Pstart,(p+1)*sizeof(long),0);			
		}
	}

	
	bsp_sync();

	bsp_pop_reg(MatrixI);

	matrix.i = MatrixI;

	//for(i=0;i<length;i++) matrix.i[i] = MatrixI[i];

	if(s==1)for(i=0;i<length;i++) printf("%d: %d=%ld\n",s,i,matrix.i[i]);

	vecfreel(MatrixI);
	//bsp_pop_reg(matrix.j);
	//bsp_pop_reg(matrix.ReValue);
	//bsp_pop_reg(matrix.Pstart);

	
	/*
	print_own_part(s,matrix);
	
	char outputname[100];
	sprintf(outputname,"output_%d-%d.txt",Options.SplitStrategy,s);

	if (!(File = fopen(outputname, "w"))) printf("Unable to open output file!\n");

	MMWriteSparseMatrix(&matrix,File,NULL,&Options);
	   	
	*/
    bsp_end();

}

int main(int argc, char **argv){
 
    bsp_init(bsp_mcl, argc, argv);
   	P = atoi(argv[1]);
    if (P>bsp_nprocs()){
        printf("Not enough processors available:");
        printf(" %d wanted, %d available\n", P, bsp_nprocs());
        exit(1);
    }
    bsp_mcl();
    exit(0);
}
