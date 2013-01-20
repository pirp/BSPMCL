#include <stdlib.h>
#include <stdio.h>
#include <Mondriaan.h>

void sparse_mult(struct sparsematrix m1, struct sparsematrix m2) {
	printf("%ld,%ld",m1.NrNzElts,m2.NrNzElts);
}


int main(int argc, char **argv){
	FILE *File;
	/* This structure will contain the testMatrix matrix. */
	struct sparsematrix testMatrix;

	/* Open the testMatrix matrix file. */
	if (!(File = fopen("test_matrix.mtx", "r")))
	{
		printf("Unable to open test_matrix!\n");
		return EXIT_FAILURE;
	}
	
	/* Read it from the file. */
	if (!MMReadSparseMatrix(File, &testMatrix))
	{
		printf("Unable to read test_matrix!\n");
		fclose(File);
		return EXIT_FAILURE;
	}
	fclose(File);
	int j=0;
	long int nonzero = testMatrix.NrNzElts;
	for(j=0;j<nonzero;j++) printf("%d: (%ld,%ld)=%1.0f\n",j,testMatrix.i[j]+1,testMatrix.j[j]+1,testMatrix.ReValue[j]);
	printf("\n");
	sparse_mult(testMatrix,testMatrix);
}


