#include "BSPedupack1.01/bspedupack.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Mondriaan.h>
#include "sp_utilities.c"

#define Niter 2

int P;
char inputname[100];

void print_matrix(int s,struct sparsematrix matrix){
	int k;
	printf("nz: %ld\n", matrix.NrNzElts);
	for(k=0;k<matrix.NrNzElts;k++){
		printf("%d: (%ld,%ld)=%f\n",s,matrix.i[k]+1,matrix.j[k]+1,matrix.ReValue[k]);
	}
}

void print_own_part(int s, struct sparsematrix matrix){
	int k;
	int start = matrix.Pstart[s];
	int end = matrix.Pstart[s+1];
	for(k=start;k<end;k++) printf("%d: (%ld,%ld)=%f\n",s,matrix.i[k]+1,matrix.j[k]+1,matrix.ReValue[k]);
}

struct sparsematrix reorder_col_incr(struct sparsematrix matrix){
	long length = matrix.NrNzElts;
	long* I = (long*) malloc(length*sizeof(long));
	long* J = (long*) malloc(length*sizeof(long));
	double* Val = (double*) malloc(length*sizeof(double));

	int k,l;

	long* tempArray = (long*) malloc(length*SZLONG);
	for(k=0;k<length;k++) tempArray[k] = matrix.j[k];

	struct sparsematrix newmatrix;

	long* indices = QSort(tempArray,length);

	for(l=0;l<length;l++){
		k = indices[length-l-1];
		I[l] = matrix.i[k];
		J[l] = matrix.j[k];
		Val[l] = matrix.ReValue[k];
		
	}

	newmatrix.i = I;
	newmatrix.j = J;
	newmatrix.ReValue = Val;
	newmatrix.NrNzElts = length;
	free(tempArray);
	return newmatrix;
}

struct sparsematrix outerproduct_even(int s, struct sparsematrix matrixA, struct sparsematrix matrixB){
	int i,j,k = 0;
	int index_A = 0;
	long int rowB,colB,rowA,colA;
	int index_B = matrixB.Pstart[s];
	struct sparsematrix finalMatrix;
	int p_end = matrixB.Pstart[s+1];
	double** val;
	val = (double**) malloc(matrixA.m*sizeof(double*));
	for(i=0;i<matrixA.m;i++) val[i] = (double*) malloc(matrixB.n*sizeof(double));
	for(i=0;i<matrixA.m;i++) for(j=0;j<matrixB.n;j++) val[i][j] =0;
	while(index_B < p_end && index_A < matrixA.NrNzElts){
		rowB = matrixB.i[index_B];
		colA = matrixA.j[index_A];

		colB = matrixB.j[index_B];
		rowA = matrixA.i[index_A];

		//printf("considering (%ld,%ld) & (%ld,%ld)\n",rowA+1,colA+1,rowB+1,colB+1);

		if(colA < rowB){
		//printf("colA<rowB!\n");
			index_A++;
			continue;
		}
		if(rowB < colA){
		//printf("rowB<colA!\n");
			index_B++;
			continue;
		}
		int index2 = index_A;
		while(matrixA.j[index2] == rowB){
			rowA = matrixA.i[index2];
			colA = matrixA.j[index2];
			val[rowA][colB] += matrixA.ReValue[index2]*matrixB.ReValue[index_B];
			//printf("(%ld,%ld)*(%ld,%ld) = (%ld,%ld)\n",rowA+1,colA+1,rowB+1,colB+1,rowA+1,colB+1);
			index2++;
		}
		index_B++;
	}

	int count = 0;

	for(i=0;i<matrixA.m;i++) for(j=0;j<matrixB.n;j++) if(val[i][j] != 0) count++;
	long* rows = (long*) malloc(count*sizeof(long));
	long* cols = (long*) malloc(count*sizeof(long));
	double* vals = (double*) malloc(count*sizeof(double));

	
	for(i=0;i<matrixA.n;i++) for(j=0;j<matrixB.n;j++){
		if(val[i][j]!=0){
			rows[k] = i;
			cols[k] = j;
			vals[k] = val[i][j];
			k++;
		}
	}
	//for(k=0;k<count;k++) printf("(%d,%d)=%f\n",rows[k]+1,cols[k]+1,vals[k]+1);

	finalMatrix.j = cols;
	finalMatrix.i = rows;
	finalMatrix.ReValue = vals;
	finalMatrix.m = matrixA.m;
	finalMatrix.n = matrixB.n;
	finalMatrix.NrNzElts = count;


	for(i=0;i<matrixB.n;i++) free(val[i]);
	free(val);
	
	return finalMatrix;
}

struct sparsematrix outerproduct_odd(int s, struct sparsematrix matrixA, struct sparsematrix matrixB){
	int i,j,k = 0;
	int index_A = 0;
	long int rowB,colB,rowA,colA;
	int p_start = matrixB.Pstart[s];
	struct sparsematrix finalMatrix;
	int index_B = matrixB.Pstart[s+1]-1;
	double** val;
	val = (double**) malloc(matrixA.m*sizeof(double*));
	for(i=0;i<matrixA.m;i++) val[i] = (double*) malloc(matrixB.n*sizeof(double));
	for(i=0;i<matrixA.m;i++) for(j=0;j<matrixB.n;j++) val[i][j] =0;

	while(index_B >= p_start && index_A < matrixA.NrNzElts){
		rowB = matrixB.i[index_B];
		colA = matrixA.j[index_A];

		colB = matrixB.j[index_B];
		rowA = matrixA.i[index_A];

		//printf("considering (%ld,%ld) & (%ld,%ld)\n",rowA+1,colA+1,rowB+1,colB+1);

		if(colA < rowB){
		//printf("colA<rowB!\n");
			index_A++;
			continue;
		}
		if(rowB < colA){
		//printf("rowB<colA!\n");
			index_B--;
			continue;
		}
		int index2 = index_A;
		while(matrixA.j[index2] == rowB){
			rowA = matrixA.i[index2];
			colA = matrixA.j[index2];
			val[rowA][colB] += matrixA.ReValue[index2]*matrixB.ReValue[index_B];
			//printf("(%ld,%ld)*(%ld,%ld) = (%ld,%ld)\n",rowA+1,colA+1,rowB+1,colB+1,rowA+1,colB+1);
			index2++;
		}
		index_B--;
	}

	int count = 0;

	for(i=0;i<matrixA.m;i++) for(j=0;j<matrixB.n;j++) if(val[i][j] != 0) count++;
	long* rows = (long*) malloc(count*sizeof(long));
	long* cols = (long*) malloc(count*sizeof(long));
	double* vals = (double*) malloc(count*sizeof(double));

	
	for(i=0;i<matrixA.n;i++) for(j=0;j<matrixB.n;j++){
		if(val[i][j]!=0){
			rows[k] = i;
			cols[k] = j;
			vals[k] = val[i][j];
			k++;
		}
	}
	//for(k=0;k<count;k++) printf("(%d,%d)=%f\n",rows[k]+1,cols[k]+1,vals[k]+1);

	finalMatrix.j = cols;
	finalMatrix.i = rows;
	finalMatrix.ReValue = vals;
	finalMatrix.m = matrixA.m;
	finalMatrix.n = matrixB.n;
	finalMatrix.NrNzElts = count;


	for(i=0;i<matrixB.n;i++) free(val[i]);
	free(val);
	
	return finalMatrix;
}

struct sparsematrix mergeLists(long* listI, long* listJ, double* listVal, int m, int n, int listSize){
	int i,j,k;
	struct sparsematrix finalMatrix;
	double** val = (double**) malloc(m*sizeof(double*));
	for(i=0;i<m;i++) val[i] = (double*) malloc(n*sizeof(double));
	for(i=0;i<m;i++) for(j=0;j<n;j++) val[i][j] =0;

	for(k=0;k<listSize;k++){
		i = listI[k];
		j = listJ[k];
		val[i][j] += listVal[k];
	}

	int count = 0;

	for(i=0;i<m;i++) for(j=0;j<n;j++) if(val[i][j] != 0) count++;
	long* rows = (long*) malloc(count*sizeof(long));
	long* cols = (long*) malloc(count*sizeof(long));
	double* vals = (double*) malloc(count*sizeof(double));

	k=0;	
	for(i=0;i<m;i++) for(j=0;j<n;j++){
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
	finalMatrix.m = m;
	finalMatrix.n = n;
	finalMatrix.NrNzElts = count;


	for(i=0;i<n;i++) free(val[i]);
	free(val);
	
	return finalMatrix;
}

struct sparsematrix sparse_mult(int s, int p, struct sparsematrix matrixA, struct sparsematrix matrixB){

	int t;

	// computes the outerproduct of the owned rows/columns
	struct sparsematrix newmatrix = (s%2) ? outerproduct_odd(s,matrixA,matrixB) : outerproduct_even(s,matrixA,matrixB);
	
	// allocation & registration of the vector that will hold the count of nonzeros from all the procs
	long* GlobNZ = vecallocl(p);
	bsp_push_reg(GlobNZ,p*SZLONG);
	bsp_sync();

	
	// every processor communicates the size of his computed part
	for(t=0;t<p;t++) bsp_put(t,&newmatrix.NrNzElts,GlobNZ,s*SZLONG,SZLONG);
	bsp_sync();
	bsp_pop_reg(GlobNZ);

	// determining the size of the list which will contain the computed part for every processor
	long finalListSize = 0;
	for(t=0;t<p;t++) finalListSize += GlobNZ[t];

	// calculating the offset in which to put the local elements
	long offset = 0;
	for(t=0;t<s;t++) offset +=GlobNZ[t];

	// lists with all the newly computed elements
	long* finalListI = vecallocl(finalListSize);
	long* finalListJ = vecallocl(finalListSize);
	double* finalListVal = vecallocd(finalListSize);
	bsp_push_reg(finalListI,finalListSize*SZLONG);
	bsp_push_reg(finalListJ,finalListSize*SZLONG);
	bsp_push_reg(finalListVal,finalListSize*SZDBL);
	bsp_sync();

	// communication
	for(t=0;t<p;t++){
		bsp_put(t,newmatrix.i,finalListI,offset*SZLONG,newmatrix.NrNzElts*SZLONG);
		bsp_put(t,newmatrix.j,finalListJ,offset*SZLONG,newmatrix.NrNzElts*SZLONG);
		bsp_put(t,newmatrix.ReValue,finalListVal,offset*SZDBL,newmatrix.NrNzElts*SZDBL);
	}

	bsp_sync();

	bsp_pop_reg(finalListI);
	bsp_pop_reg(finalListJ);
	bsp_pop_reg(finalListVal);

	// merging duplicates in this list, returning a new sparse matrix
	struct sparsematrix matrix2 = mergeLists(finalListI,finalListJ,finalListVal,newmatrix.m,newmatrix.n,finalListSize);
	
	vecfreel(GlobNZ);
	vecfreel(finalListI);
	vecfreel(finalListJ);
	vecfreed(finalListVal);	

	return matrix2;
}

void read_from_file(char* input, struct sparsematrix *matrixA, struct sparsematrix *matrixB){
	// every processor reads matrices from file, such that arrays are allocated
    FILE *File;
    

	if (!(File = fopen(input, "r"))) printf("Unable to open the matrix!\n");
	if (!MMReadSparseMatrix(File, matrixA)) printf("Unable to read into A!\n");
	fclose(File);

	if (!(File = fopen(input, "r"))) printf("Unable to open the matrix!\n");
	if (!MMReadSparseMatrix(File, matrixB)) printf("Unable to read into B!\n");
	fclose(File);
}

void split_matrices(int s, int p, struct sparsematrix *matrixB){
	int t;
	// length of the arrays in the sparsematrix struct
	long length = matrixB->NrNzElts;

	// temporary arrays for communication
	long *MatrixBI, *MatrixBJ, *MatrixBPstart;
	double  *MatrixBReValue;
	MatrixBI = vecallocl(length);
	MatrixBJ = vecallocl(length);
	MatrixBReValue = vecallocd(length);
	MatrixBPstart = vecallocl(p+1);
	
	//allocation of Pstart for the processors that will not get it from mondriaan
	if (s!=0) matrixB->Pstart = vecallocl(p+1);

	// registration of the temp arrays
	bsp_push_reg(MatrixBI,length*SZLONG);
	bsp_push_reg(MatrixBJ,length*SZLONG);
	bsp_push_reg(MatrixBReValue,length*SZDBL);
	bsp_push_reg(MatrixBPstart,length*SZLONG);

	bsp_sync();


	if(s==0){
		// only proc 0 uses mondriaan and splits the matrix 

		//print_matrix(s,*matrixB);
		struct opts OptionsB;
		if (!SetOptionsFromFile(&OptionsB, "Mondriaan.defaults")) printf("Unable to set options from disk!\n");
		OptionsB.SplitStrategy = 4;
		if (!ApplyOptions(&OptionsB)) printf("Invalid options!\n");
		if (!DistributeMatrixMondriaan(matrixB, p, 0.03, &OptionsB, NULL)) printf("Unable to distribute!\n");
		
		// communication
		for(t=0;t<p;t++){
			bsp_put(t,matrixB->i,MatrixBI,0,length*SZLONG);			
			bsp_put(t,matrixB->j,MatrixBJ,0,length*SZLONG);			
			bsp_put(t,matrixB->ReValue,MatrixBReValue,0,length*SZDBL);			
			bsp_put(t,matrixB->Pstart,MatrixBPstart,0,(p+1)*SZLONG);
		}
	}

	
	bsp_sync();

	// B is updated with the new information from p0

	for(t=0;t<length;t++){
		matrixB->i[t] = MatrixBI[t];
		matrixB->j[t] = MatrixBJ[t];
		matrixB->ReValue[t] = MatrixBReValue[t];
	}

	for(t=0;t<p+1;t++) matrixB->Pstart[t] = MatrixBPstart[t];

	// registration of the temp arrays not needed anymore
	bsp_pop_reg(MatrixBI);
	bsp_pop_reg(MatrixBJ);
	bsp_pop_reg(MatrixBReValue);
	bsp_pop_reg(MatrixBPstart);

	// temp arrays not required anymore are freed from memory
	vecfreel(MatrixBI);
	vecfreel(MatrixBJ);
	vecfreed(MatrixBReValue);
	vecfreel(MatrixBPstart);
}

struct sparsematrix iteration(int s, int p, struct sparsematrix matrixA, struct sparsematrix matrixB){
	split_matrices(s,p,&matrixB);
	struct sparsematrix matrix = sparse_mult(s,p,matrixA,matrixB);
	matrix = inflateMatrix(matrix,2);
	return matrix;
}


void bsp_mcl(){
	int p,s,t;
	int i,j;

	bsp_begin(P);

	p= bsp_nprocs();
    s= bsp_pid();

    struct sparsematrix matrixA;
    struct sparsematrix matrixB;
    struct sparsematrix matrix;

    read_from_file(inputname,&matrixA,&matrixB);


    for(i=0;i<Niter;i++){
    	matrixA = normalize_rows(matrixA);
    	matrixB = normalize_rows(matrixB);
    	if(s==0) printf("iteration: %d - %ld\n",i,matrixB.NrNzElts);
    	matrix = iteration(s,p,matrixA,matrixB);

  
	    matrixB.i = matrix.i;
	    matrixB.j = matrix.j;
	    matrixB.ReValue = matrix.ReValue;
	    matrixB.NrNzElts = matrix.NrNzElts;

	    matrixA = reorder_col_incr(matrix);
    }

  	//if(s==0) print_matrix(s,matrix);

    bsp_end();
}

int main(int argc, char **argv){
 
    bsp_init(bsp_mcl, argc, argv);
   	P = atoi(argv[1]);
   	sprintf(inputname,"%s.mtx",argv[2]);
    if (P>bsp_nprocs()){
        printf("Not enough processors available:");
        printf(" %d wanted, %d available\n", P, bsp_nprocs());
        exit(1);
    }
    bsp_mcl();
    exit(0);
}
