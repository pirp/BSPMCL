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


	for(k=0;k<matrix.NrNzElts;k++){
		col = matrix.j[k];
		value = matrix.ReValue[k];
		matrix.ReValue[k] = value/sumValues[col];
	}

	return matrix;

}

struct sparsematrix normalize_rows(struct sparsematrix matrix){
	int numberRows = matrix.m;
	double* sumValues = (double*) malloc(numberRows*sizeof(double));
	int row,k;
	double value;

	for(k=0;k<numberRows;k++) sumValues[k] = 0;

	for(k=0;k<matrix.NrNzElts;k++){
		row = matrix.i[k];
		value = matrix.ReValue[k];
		sumValues[row] += value;
	}


	for(k=0;k<matrix.NrNzElts;k++){
		row = matrix.i[k];
		value = matrix.ReValue[k];
		matrix.ReValue[k] = value/sumValues[row];
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
}