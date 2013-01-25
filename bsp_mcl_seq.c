#include <stdlib.h>
#include <stdio.h>
#include <Mondriaan.h>

void sparse_mult(struct sparsematrix m1, struct sparsematrix m2) {
	int i,j,k,r,s;
	int row, col,row2,col2;
	double** val;
	val = (double**) malloc(m1.n*sizeof(double*));
	for(i=0;i<m1.n;i++) val[i] = (double*) malloc(m1.n*sizeof(double));
	
	for(i=0;i<m1.n;i++) for(j=0;j<m1.n;j++) val[i][j] =0;
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
	/*
	for(i=0;i<m1.n;i++){
		for(j=0;j<m1.n;j++) printf("%1.0f ",val[i][j]);
		printf("\n");
	}
	*/
	printf("\n\n\n");

	int count =0;
	for(i=0;i<m1.n;i++) for(j=0;j<m1.n;j++) if(val[i][j] != 0) count++;
	//printf("nonzeros=%d\n\n",count);
	

	//int* rows,cols,vals;
	int* rows = (int*) malloc(count*sizeof(int));
	int* cols = (int*) malloc(count*sizeof(int));
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
	/*
	for(k=0;k<count;k++){
		printf("(%d,%d)=%f\n",rows[k]+1,cols[k]+1,vals[k]);
	}
*/
	for(i=0;i<m1.n;i++) free(val[i]);
	free(val);
}


int main(int argc, char **argv){
	FILE* File;
	struct sparsematrix testMatrix;

	if (!(File = fopen("bcsstk08.mtx", "r")))
	{
		printf("Unable to open test_matrix!\n");
		return EXIT_FAILURE;
	}
	
	if (!MMReadSparseMatrix(File, &testMatrix))
	{
		printf("Unable to read test_matrix!\n");
		fclose(File);
		return EXIT_FAILURE;
	}

	sparse_mult(testMatrix,testMatrix);
	
	fclose(File);
}


