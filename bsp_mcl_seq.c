#include <stdlib.h>
#include <stdio.h>
#include <Mondriaan.h>
#include <math.h>
#include "sp_utilities.c"
#include <time.h>

#define Niter 15

void print_matrix(struct sparsematrix matrix){
	int k;
	for(k=0;k<matrix.NrNzElts;k++){
		printf("(%ld,%ld)=%f\n", matrix.i[k],matrix.j[k],matrix.ReValue[k]);
	}
}

struct sparsematrix sparse_mult(struct sparsematrix m1, struct sparsematrix m2) {
	int i,j,k,r,s;
	struct sparsematrix finalMatrix;
	int row, col,row2,col2;
	double** val;
	val = (double**) malloc(m1.n*sizeof(double*));
	for(i=0;i<m1.n;i++) val[i] = (double*) malloc(m1.n*sizeof(double));
	
	for(i=0;i<m1.n;i++) for(j=0;j<m2.n;j++) val[i][j] =0;
	for(k=0;k<m1.NrNzElts;k++){
		row = m1.i[k];
		col = m1.j[k];
		for(r=0;r<m2.NrNzElts;r++){
			row2 = m2.i[r];
			col2 = m2.j[r];
			if(row2 == col){
				val[row][col2] += m1.ReValue[k]*m2.ReValue[r];
			}
		}
	}

	//for(i=0;i<m1.n;i++) for(j=0;j<m1.n;j++) if(val[i][j] < 0.001) val[i][j] =0;

	int count =0;
	for(i=0;i<m1.n;i++) for(j=0;j<m1.n;j++) if(val[i][j] != 0) count++;
	
	long* rows = (long*) malloc(count*sizeof(long));
	long* cols = (long*) malloc(count*sizeof(long));
	double* vals = (double*) malloc(count*sizeof(double));

	k=0;
	for(i=0;i<m1.n;i++) for(j=0;j<m1.n;j++){
		if(val[i][j]!=0){
			rows[k] = i;
			cols[k] = j;
			vals[k] = val[i][j];
			k++;
		}
	}

	finalMatrix.j = cols;
	finalMatrix.i = rows;
	finalMatrix.ReValue = vals;
	finalMatrix.m = m1.m;
	finalMatrix.n = m2.n;
	finalMatrix.NrNzElts = count;


	for(i=0;i<m1.n;i++) free(val[i]);
	free(val);
	
	return finalMatrix;
}

struct sparsematrix power_matrix(struct sparsematrix matrix, int e){
	int i=1;
	while(i<e){
		matrix = sparse_mult(matrix,matrix);
		i++;
	}
	return matrix;
}

struct sparsematrix iteration(struct sparsematrix matrix, int r, int e){
	matrix = normalize_rows(matrix);
	matrix = power_matrix(matrix,e);
	matrix = inflateMatrix(matrix,r);
	return matrix;
}

int main(int argc, char **argv){
	FILE* File;
	struct sparsematrix matrix;
	char inputname[100];
	sprintf(inputname,"%s",argv[1]);

	if (!(File = fopen(inputname, "r"))) printf("Unable to open input matrix! Make sure it's the first parameter to the program\n");
	if (!MMReadSparseMatrix(File, &matrix))	printf("Unable to read input matrix!\n");
	fclose(File);

	int k;
	clock_t time_start,time0,time1;
	time_start = clock();
	time0 = time_start;

	//printf("iteration: %d \t %ld \n",k+1,matrix.NrNzElts);
	
	for(k=0;k<Niter;k++){
		matrix = iteration(matrix,2,2);
		time1 = clock();
		//printf("iteration: %d \t %ld \t %f \n",k+1,matrix.NrNzElts,(double) (time1-time0)/CLOCKS_PER_SEC);
		time0=time1;
		if (matrix.NrNzElts ==0) break;
	}
	printf("%f \n",(double) (time1-time_start)/CLOCKS_PER_SEC);
	

	char outputname[100];
	sprintf(outputname,"%s.txt",argv[1]);

	if (!(File = fopen(outputname, "w")))
	{
		printf("Unable to open output file!\n");
		return EXIT_FAILURE;
	}

	for(k=0;k<matrix.NrNzElts;k++) fprintf(File,"%ld %ld %f\n",matrix.i[k]+1,matrix.j[k]+1,matrix.ReValue[k]);
	fclose(File);
	
}


