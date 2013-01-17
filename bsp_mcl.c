/*
This is a small example program which illustrates interfacing with the Mondriaan library.

It can be compiled using:

gcc bsp_mcl.c -IMondriaan3.11/src/include -LMondriaan3.11/src/lib -lMondriaan3 -lm

*/

#include <stdlib.h>
#include <stdio.h>
/* Make sure to include the Mondriaan headers. */
#include <Mondriaan.h>

int main(int argc, char **argv)
{
	/* This will contain the Mondriaan options. */
	struct opts Options;
	/* This file pointer will be the opened Matrix Market file. */
	FILE *File;
	/* This structure will contain the testMatrix matrix. */
	struct sparsematrix testMatrix;
	/* Variables used for calculating the communication volume. */
	long ComVolumeRow, ComVolumeCol, Dummy;
	/* Variables used to calculate the imbalance. */
	long MaxNrNz;
	double Imbalance;
	
	/* Set the default options. */
	SetDefaultOptions(&Options);
	
	/* We can also read options from disk. */
	/* if (!SetOptionsFromFile(&Options, "options.txt")) printf("Unable to set options from disk!\n"); */
	
	/* If we are done setting the options, we check and apply them. */
	if (!ApplyOptions(&Options))
	{
		printf("Invalid options!\n");
		return EXIT_FAILURE;
	}
	
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
	
	/* Distribute the matrix over two processors with an allowed imbalance of 3% and the options provided above. */
	if (!DistributeMatrixMondriaan(&testMatrix, 4, 0.03, &Options, NULL))
	{
		printf("Unable to distribute test_matrix!\n");
		return EXIT_FAILURE;
	}
	
	/* Calculate the communication volume. */
	CalcCom(&testMatrix, NULL, ROW, &ComVolumeRow, &Dummy, &Dummy, &Dummy, &Dummy);
	CalcCom(&testMatrix, NULL, COL, &ComVolumeCol, &Dummy, &Dummy, &Dummy, &Dummy);
	
	/* Calculate the imbalance, making use of the fact that we only distributed the nonzeros over two processors. */
	MaxNrNz = MAX(testMatrix.Pstart[2] - testMatrix.Pstart[1], testMatrix.Pstart[1] - testMatrix.Pstart[0]);
	Imbalance = (double)(2*MaxNrNz - testMatrix.NrNzElts)/(double)testMatrix.NrNzElts;
	
	if (Imbalance > 0.03)
	{
		printf("Imbalance is too large!\n");
		return EXIT_FAILURE;
	}
	
	/* Display information about this partitioning. */
	printf("Succesfully distributed %ld nonzeros over two processors: %ld are assigned to processor 0 and %ld to processor 1.\n", testMatrix.NrNzElts, testMatrix.Pstart[1] - testMatrix.Pstart[0], testMatrix.Pstart[2] - testMatrix.Pstart[1]);
	printf("This distribution has a total communication volume equal to %ld and imbalance equal to %.1f%%.\n", ComVolumeRow + ComVolumeCol, 100.0*Imbalance);
	
	char output[MAX_WORD_LENGTH]; /* filename of the output */

	sprintf(output, "test-P%d", testMatrix.NrProcs);
	File = fopen(output, "w");
	MMWriteSparseMatrix(&testMatrix, File, NULL, &Options);
	fclose(File);

	/* Free matrix data. */
	MMDeleteSparseMatrix(&testMatrix);
	
	return EXIT_SUCCESS;
}

