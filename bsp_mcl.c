#include "BSPedupack1.01/bspedupack.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Mondriaan.h>

int P;
char inputname[100];

void print_matrix(int s,struct sparsematrix matrix){
	int k;
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

void bsp_mcl(){
	int p,s,t;
	int i,j;

	bsp_begin(P);

	p= bsp_nprocs();
    s= bsp_pid();

    // every processor reads the matrix from file, such that arrays are allocated
    FILE *File;
    struct sparsematrix matrixA;
    struct sparsematrix matrixB;

	if (!(File = fopen(inputname, "r"))) printf("Unable to open the matrix!\n");
	if (!MMReadSparseMatrix(File, &matrixA)) printf("Unable to read into A!\n");
	fclose(File);

	if (!(File = fopen(inputname, "r"))) printf("Unable to open the matrix!\n");
	if (!MMReadSparseMatrix(File, &matrixB)) printf("Unable to read into B!\n");
	fclose(File);

	// length of the arrays in the sparsematrix struct
	long length = matrixA.NrNzElts;

	// temporary arrays for communication
	//long *MatrixAI, *MatrixAJ, *MatrixAPstart;
	long *MatrixBI, *MatrixBJ, *MatrixBPstart;
	//double  *MatrixAReValue;
	double  *MatrixBReValue;
	/*MatrixAI = vecallocl(length);
	MatrixAJ = vecallocl(length);
	MatrixAReValue = vecallocd(length);
	MatrixAPstart = vecallocl(p+1);*/
	MatrixBI = vecallocl(length);
	MatrixBJ = vecallocl(length);
	MatrixBReValue = vecallocd(length);
	MatrixBPstart = vecallocl(p+1);
	if (s!=0) {
		//matrixA.Pstart = vecallocl(p+1);
		matrixB.Pstart = vecallocl(p+1);
	}

	// registration of the temp arrays
	/*bsp_push_reg(MatrixAI,length*SZLONG);
	bsp_push_reg(MatrixAJ,length*SZLONG);
	bsp_push_reg(MatrixAReValue,length*SZDBL);
	bsp_push_reg(MatrixAPstart,length*SZLONG);*/
	bsp_push_reg(MatrixBI,length*SZLONG);
	bsp_push_reg(MatrixBJ,length*SZLONG);
	bsp_push_reg(MatrixBReValue,length*SZDBL);
	bsp_push_reg(MatrixBPstart,length*SZLONG);

	
	bsp_sync();

	if(s==0){
		// only proc 0 uses mondriaan and splits the matrix 
		//struct opts OptionsA;
		struct opts OptionsB;

		//if (!SetOptionsFromFile(&OptionsA, "Mondriaan.defaults")) printf("Unable to set options from disk!\n");
		if (!SetOptionsFromFile(&OptionsB, "Mondriaan.defaults")) printf("Unable to set options from disk!\n");

		//OptionsA.SplitStrategy = 5;
		OptionsB.SplitStrategy = 4;
		//if (!ApplyOptions(&OptionsA)) printf("Invalid options!\n");
		if (!ApplyOptions(&OptionsB)) printf("Invalid options!\n");


		//if (!DistributeMatrixMondriaan(&matrixA, p, 0.03, &OptionsA, NULL)) printf("Unable to distribute!\n");
		if (!DistributeMatrixMondriaan(&matrixB, p, 0.03, &OptionsB, NULL)) printf("Unable to distribute!\n");

		
		// communication
		for(t=0;t<p;t++){
			/*bsp_put(t,matrixA.i,MatrixAI,0,length*SZLONG);			
			bsp_put(t,matrixA.j,MatrixAJ,0,length*SZLONG);			
			bsp_put(t,matrixA.ReValue,MatrixAReValue,0,length*SZDBL);			
			bsp_put(t,matrixA.Pstart,MatrixAPstart,0,(p+1)*SZLONG);*/
			bsp_put(t,matrixB.i,MatrixBI,0,length*SZLONG);			
			bsp_put(t,matrixB.j,MatrixBJ,0,length*SZLONG);			
			bsp_put(t,matrixB.ReValue,MatrixBReValue,0,length*SZDBL);			
			bsp_put(t,matrixB.Pstart,MatrixBPstart,0,(p+1)*SZLONG);
		}
		/* char outputname[100];
		sprintf(outputname,"outputB_%d-%d.txt",OptionsB.SplitStrategy,s);
		if (!(File = fopen(outputname, "w"))) printf("Unable to open output file!\n");
		MMWriteSparseMatrix(&matrixB,File,NULL,&OptionsB);
		fclose(File);*/
	}

	
	bsp_sync();

	// sparsematrix is updated with the new information from p0

	for(t=0;t<length;t++){
		/*matrixA.i[t] = MatrixAI[t];
		matrixA.j[t] = MatrixAJ[t];
		matrixA.ReValue[t] = MatrixAReValue[t];*/
		matrixB.i[t] = MatrixBI[t];
		matrixB.j[t] = MatrixBJ[t];
		matrixB.ReValue[t] = MatrixBReValue[t];
	}


	for(t=0;t<p+1;t++){	
		//matrixA.Pstart[t] = MatrixAPstart[t];
		matrixB.Pstart[t] = MatrixBPstart[t];
	}

	//if (s==0) for(t=0;t<p+1;t++) printf("%d: %d %ld\n", s, t, matrixA.Pstart[t]);

	// registration of the temp arrays not needed anymore
	/*bsp_pop_reg(MatrixAI);
	bsp_pop_reg(MatrixAJ);
	bsp_pop_reg(MatrixAReValue);
	bsp_pop_reg(MatrixAPstart);*/
	bsp_pop_reg(MatrixBI);
	bsp_pop_reg(MatrixBJ);
	bsp_pop_reg(MatrixBReValue);
	bsp_pop_reg(MatrixBPstart);

	// temp arrays are freed from memory
	/*vecfreel(MatrixAI);
	vecfreel(MatrixAJ);
	vecfreed(MatrixAReValue);
	vecfreel(MatrixAPstart);*/
	vecfreel(MatrixBI);
	vecfreel(MatrixBJ);
	vecfreed(MatrixBReValue);
	vecfreel(MatrixBPstart);

	struct sparsematrix newmatrix = (s%2) ? outerproduct_odd(s,matrixA,matrixB) : outerproduct_even(s,matrixA,matrixB);
	
	long* GlobNZ = vecallocl(p); bsp_push_reg(GlobNZ,p*SZLONG);
	bsp_sync();
	
	for(t=0;t<p;t++) bsp_put(t,&newmatrix.NrNzElts,GlobNZ,s*SZLONG,SZLONG);

	bsp_sync();

	bsp_pop_reg(GlobNZ);

	long finalListSize = 0;
	for(t=0;t<p;t++) finalListSize += GlobNZ[t];
	long offset = 0;
	for(t=0;t<s;t++) offset +=GlobNZ[t];

	long* finalListI = vecallocl(finalListSize);
	long* finalListJ = vecallocl(finalListSize);
	double* finalListVal = vecallocd(finalListSize);


	bsp_push_reg(finalListI,finalListSize*SZLONG);
	bsp_push_reg(finalListJ,finalListSize*SZLONG);
	bsp_push_reg(finalListVal,finalListSize*SZDBL);
	bsp_sync();


	for(t=0;t<p;t++){
		bsp_put(t,newmatrix.i,finalListI,offset*SZLONG,newmatrix.NrNzElts*SZLONG);
		bsp_put(t,newmatrix.j,finalListJ,offset*SZLONG,newmatrix.NrNzElts*SZLONG);
		bsp_put(t,newmatrix.ReValue,finalListVal,offset*SZDBL,newmatrix.NrNzElts*SZDBL);
	}

	bsp_sync();

	struct sparsematrix matrix2 = mergeLists(finalListI,finalListJ,finalListVal,newmatrix.m,newmatrix.n,finalListSize);

	if(s==0) print_matrix(s,matrix2);



	bsp_pop_reg(finalListI);
	bsp_pop_reg(finalListJ);
	bsp_pop_reg(finalListVal);
	
	vecfreel(GlobNZ);
	vecfreel(finalListI);
	vecfreel(finalListJ);
	vecfreed(finalListVal);

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
