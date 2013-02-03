#define DIV 0
#define MOD 1
#define NITERS 10
#define STRLEN 100
typedef struct {int i,j;} indexpair;


void bspinput2triple(int p, int s, int *pnA, int *pnz, 
                     int **pia, int **pja, double **pa){
  
    /* This function reads a sparse matrix in distributed
       Matrix Market format without the banner line
       from the input file and distributes
       matrix triples to the processors.
       The input consists of one line
           m n nz p  (number of rows, columns, nonzeros, processors)
       followed by p+1 lines with the starting numbers
       of the processor parts
           Pstart[0]
           Pstart[1]
           ...
           Pstart[p]
       which means that processor q will get all nonzeros
       numbered Pstart[q]..Pstart[q+1]-1.
       This is followed by nz lines in the format
           i j a     (row index, column index, numerical value).
       The input indices are assumed by Matrix Market to start
       counting at one, but they are converted to start from zero.
       The triples are stored into three arrays ia, ja, a,
       in arbitrary order.
       
       Input:
       p is the number of processors.
       s is the processor number, 0 <= s < p.

       Output:
       nA is the global matrix size.
       nz is the number of local nonzeros.
       a[k] is the numerical value of the k'th local nonzero,
            0 <= k < nz.
       ia[k] is the global row index of the  k'th local nonzero.
       ja[k] is the global column index.
    */

    int pA, mA, nA, nzA, nz, q, nzq, k, tagsz, status, *Pstart, *ia, *ja;
    double value, *a;
    indexpair t;
    char filename[STRLEN];
    FILE *fp;
    
    Pstart= vecalloci(p+1);
    bsp_push_reg(&nA,SZINT);
    bsp_push_reg(&nz,SZINT);
    tagsz= sizeof(indexpair);
    bsp_set_tagsize(&tagsz);
    bsp_sync();

    if (s==0){
        /* Open the matrix file and read the header */
        printf("Please enter the filename of the matrix distribution\n");
        scanf("%s",filename);
        fp=fopen(filename,"r");

        /* A is an mA by nA matrix with nzA nonzeros
           distributed over pA processors. */
        fscanf(fp,"%d %d %d %d\n", &mA, &nA, &nzA, &pA);
        if(pA!=p)
            bsp_abort("Error: p not equal to p(A)\n"); 
        if(mA!=nA)
            bsp_abort("Error: matrix is not square");

        for (q=0; q<=p; q++)
            fscanf(fp,"%d\n", &Pstart[q]);
        for (q=0; q<p; q++){
            bsp_put(q,&nA,&nA,0,SZINT);
            nzq= Pstart[q+1]-Pstart[q];
            bsp_put(q,&nzq,&nz,0,SZINT);
        }
    }
    bsp_sync();


    /* Handle the processors one at a time.
       This saves buffer memory, at the expense of p-1 extra syncs.
       Buffer memory needed for communication is at most the maximum
       amount of memory a processor needs to store its vector components. */

    a= vecallocd(nz+1);
    ia= vecalloci(nz+1);  
    ja= vecalloci(nz+1);

    for (q=0; q<p; q++){      
        if (s==0){
            /* Read the nonzeros from the matrix file and
               send them to their destination */
            for (k=Pstart[q]; k<Pstart[q+1]; k++){
                fscanf(fp,"%d %d %lf\n", &t.i, &t.j, &value);
                /* Convert indices to range 0..n-1, 
                   assuming it was 1..n */
                t.i--;
                t.j--;
                /* Send a triple to P(q). Tag is a pair (i,j).
                   Payload is a numerical value */
                bsp_send(q,&t,&value,SZDBL);
            }
        }
        bsp_sync();
        
        if (s==q){
            /* Store the received nonzeros */
            for(k=0; k<nz; k++){
                bsp_get_tag(&status,&t);
                ia[k]= t.i;
                ja[k]= t.j;
                bsp_move(&a[k],SZDBL);
            }
        }
    }

    *pnA= nA;
    *pnz= nz;
    *pa= a;
    *pia= ia;
    *pja= ja;
    if (s==0)
        fclose(fp);
    bsp_pop_reg(&nz);
    bsp_pop_reg(&nA);
    vecfreei(Pstart);
    bsp_sync();
    
} /* end bspinput2triple */

int key(int i, int radix, int keytype){
   /* This function computes the key of an index i
      according to the keytype */
   
       if (keytype==DIV)
           return i/radix;
       else /* keytype=MOD */
           return i%radix;
           
} /* end key */

void sort(int n, int nz, int *ia, int *ja, double *a,int radix, int keytype){
   /* This function sorts the nonzero elements of an n by n sparse
      matrix A stored in triple format in arrays ia, ja, a.
      The sort is by counting. 
      If keytype=DIV, the triples are sorted by increasing value of
      ia[k] div radix.
      if keytype=MOD, the triples are sorted by increasing value of
      ia[k] mod radix.
      The sorting is stable: ties are decided so that the original
      precedences are maintained. For a complete sort by increasing
      index ia[k], this function should be called twice:
      first with keytype=MOD, then with keytype=DIV.
      
      Input: 
      n is the global size of the matrix.
      nz is the local number of nonzeros.
      a[k] is the numerical value of the k'th nonzero of the
           sparse matrix A, 0 <= k < nz.
      ia[k] is the global row index of the k'th nonzero.
      ja[k] is the global column index of the k'th nonzero.
      radix >= 1.
      
      Output: ia, ja, a in sorted order.
   */
   
   int key(int i, int radix, int keytype);

   int *ia1, *ja1, nbins, *startbin, *lengthbin, r, k, newk;
   double *a1;
   
   ia1= vecalloci(nz);
   ja1= vecalloci(nz); 
   a1 = vecallocd(nz);
   
   /* Allocate bins */
   if (keytype==DIV)
       nbins= (n%radix==0 ? n/radix : n/radix+1);
   else if (keytype==MOD)
       nbins= radix;
   startbin= vecalloci(nbins);
   lengthbin= vecalloci(nbins);
       
   /* Count the elements in each bin */
   for (r=0; r<nbins; r++)
       lengthbin[r]= 0;
   for (k=0; k<nz; k++){
       r= key(ia[k],radix,keytype);
       lengthbin[r]++;
   }
    
   /* Compute the starting positions */
   startbin[0]= 0;
   for (r=1; r<nbins; r++)
       startbin[r]= startbin[r-1] + lengthbin[r-1];
       
   /* Enter the elements into the bins in temporary arrays (ia1,ja1,a1) */
   for (k=0; k<nz; k++){
       r= key(ia[k],radix,keytype);
       newk= startbin[r];
       ia1[newk]= ia[k];
       ja1[newk]= ja[k];
       a1[newk] = a[k];
       startbin[r]++;
   }
  
   /* Copy the elements back to the orginal arrays */
   for (k=0; k<nz; k++){
       ia[k]= ia1[k];
       ja[k]= ja1[k];
       a[k] = a1[k];
   }
   
   vecfreei(lengthbin);
   vecfreei(startbin);
   vecfreed(a1);
   vecfreei(ja1);
   vecfreei(ia1);
   
} /* end sort */

void triple2icrs(int n, int nz, int *ia,  int *ja, double *a,int *pnrows, int *pncols,int **prowindex, int **pcolindex){
    /* This function converts a sparse matrix A given in triple
       format with global indices into a sparse matrix in
       incremental compressed row storage (ICRS) format with 
       local indices.

       The conversion needs time and memory O(nz + sqrt(n))
       on each processor, which is O(nz(A)/p + n/p + p).
       
       Input:
       n is the global size of the matrix.
       nz is the local number of nonzeros.
       a[k] is the numerical value of the k'th nonzero
            of the sparse matrix A, 0 <= k <nz.
       ia[k] is the global row index of the k'th nonzero.
       ja[k] is the global column index of the k'th nonzero.
  
       Output:
       nrows is the number of local nonempty rows
       ncols is the number of local nonempty columns
       rowindex[i] is the global row index of the i'th
                   local row, 0 <= i < nrows.
       colindex[j] is the global column index of the j'th
                   local column, 0 <= j < ncols.
       a[k] is the numerical value of the k'th local nonzero of the
            sparse matrix A, 0 <= k < nz. The array is sorted by
            row index, ties being decided by column index.
       ia[k] = inc[k] is the increment in the local column index of the
              k'th local nonzero, compared to the column index of the
              (k-1)th nonzero, if this nonzero is in the same row;
              otherwise, ncols is added to the difference.
              By convention, the column index of the -1'th nonzero is 0.
   */
    
   void sort(int n, int nz, int *ia, int *ja, double *a,
             int radix, int keytype);

   int radix, i, iglob, iglob_last, j, jglob, jglob_last, k, inck,
       nrows, ncols, *rowindex, *colindex;
   
   /* radix is the smallest power of two >= sqrt(n)
      The div and mod operations are cheap for powers of two.
      A radix of about sqrt(n) minimizes memory and time. */

   for (radix=1; radix*radix<n; radix *= 2)
       ;
   
   /* Sort nonzeros by column index */
   sort(n,nz,ja,ia,a,radix,MOD);
   sort(n,nz,ja,ia,a,radix,DIV);
   
   /* Count the number of local columns */
   ncols= 0;
   jglob_last= -1;
   for(k=0; k<nz; k++){
       jglob= ja[k];
       if(jglob!=jglob_last)
           /* new column index */
           ncols++;
       jglob_last= jglob;
   }
   colindex= vecalloci(ncols);
   
   /* Convert global column indices to local ones.
      Initialize colindex */
   j= 0;
   jglob_last= -1;
   for(k=0; k<nz; k++){
       jglob= ja[k];
       if(jglob!=jglob_last){
           colindex[j]= jglob;
           j++;
       }
       ja[k]= j-1; /* local index of last registered column */
       jglob_last= jglob;
   }
   
   /* Sort nonzeros by row index using radix-sort */
   sort(n,nz,ia,ja,a,radix,MOD);
   sort(n,nz,ia,ja,a,radix,DIV);

   /* Count the number of local rows */
   nrows= 0;
   iglob_last= -1;
   for(k=0; k<nz; k++){
       iglob= ia[k];
       if(iglob!=iglob_last)
           /* new row index */
           nrows++;
       iglob_last= iglob;
   }
   rowindex= vecalloci(nrows);
                              
   /* Convert global row indices to local ones.
      Initialize rowindex and inc */
   i= 0;
   iglob_last= -1;
   for(k=0; k<nz; k++){
       if (k==0)
           inck= ja[k];
       else
           inck= ja[k] - ja[k-1];
       iglob= ia[k]; 
       if(iglob!=iglob_last){
           rowindex[i]= iglob;
           i++;
           if(k>0)
               inck += ncols;
       } 
       ia[k]= inck; /* ia is used to store inc */
       iglob_last= iglob;
   }
   if (nz==0)
       ia[nz]= 0;
   else 
       ia[nz]= ncols - ja[nz-1];
   ja[nz]= 0;                                                  
   a[nz]= 0.0;     
   
   *pncols= ncols;
   *pnrows= nrows;
   *prowindex= rowindex;
   *pcolindex= colindex;
   
} /* end triple2icrs */
