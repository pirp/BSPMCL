#include <stdlib.h>
#include <stdio.h>
#include <Mondriaan.h>
#include <math.h>

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
	/*
	for(k=0;k<count;k++){
		printf("(%d,%d)=%f\n",rows[k]+1,cols[k]+1,vals[k]);
	}
*/	
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

struct sparsematrix normalize_columns(struct sparsematrix matrix){
	int numberCols = matrix.n;
	double* sumValues = (double*) malloc(numberCols*sizeof(double));
	int col,k;
	double value;

	for(k=0;k<numberCols;k++) sumValues[k] = 0;

	for(k=0;k<matrix.NrNzElts;k++){
		col = matrix.j[k];
		value = matrix.ReValue[k];
		sumValues[col] += value;
	}
/*
	printf("highest values:\n");
	for(k=0;k<numberCols;k++) printf("%d: %f\n",k,sumValues[k]);
*/

	for(k=0;k<matrix.NrNzElts;k++){
		col = matrix.j[k];
		value = matrix.ReValue[k];
		matrix.ReValue[k] = value/sumValues[col];
		//printf("(%d,%d)=%f\n",matrix.i[k]+1,matrix.j[k]+1,matrix.ReValue[k]);
	}

	return matrix;

}

struct sparsematrix inflateMatrix(struct sparsematrix matrix, int r){
	int k;
	double value;
	for(k=0;k<matrix.NrNzElts;k++){
		value = matrix.ReValue[k];	
		matrix.ReValue[k] = pow(value,r);
	}
	return matrix;
};

int main(int argc, char **argv){
	FILE* File;
	struct sparsematrix testMatrix;

	if (!(File = fopen("test_matrix.mtx", "r")))
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
	struct sparsematrix matrix;

	int k;
	//newMatrix = sparse_mult(testMatrix,testMatrix);
	matrix = normalize_columns(testMatrix);
	matrix = inflateMatrix(testMatrix,2);

	//for(k=0;k<matrix.NrNzElts;k++) printf("(%d,%d)=%f\n",matrix.i[k]+1,matrix.j[k]+1,matrix.ReValue[k]);
	
	fclose(File);
}


