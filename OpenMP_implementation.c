/*
 bitonic.c 

 This file contains two different implementations of the bitonic sort
        recursive  version :  recBitonicSort()
        imperative version :  impBitonicSort() 
 

 The bitonic sort is also known as Batcher Sort. 
 For a reference of the algorithm, see the article titled 
 Sorting networks and their applications by K. E. Batcher in 1968 


 The following codes take references to the codes avaiable at 

 http://www.cag.lcs.mit.edu/streamit/results/bitonic/code/c/bitonic.c

 http://www.tools-of-computing.com/tc/CS/Sorts/bitonic_sort.htm

 http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/bitonicen.htm 
 */

/* 
------- ---------------------- 
   Nikos Pitsianis, Duke CS 
-----------------------------
*/

/* 
-------------------------------
        STUDENT:
    PIPERIDIS ANESTIS
       AEM:8689
-------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>


struct timeval startwtime, endwtime;
double seq_time;


int N;          // data array size
int *a;         // data array to be sorted

const int ASCENDING  = 1;
const int DESCENDING = 0;

int NUM_OF_THREADS;//2^p threads

void init(void);
void print(void);
void sort(void);
void test(void);
inline void exchange(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSort(int lo, int cnt, int dir);
void impBitonicSort(void);

//function useful fort the use of quicksort
int cmpfunc_asc (const void * a, const void * b);
int cmpfunc_dsc (const void * a, const void * b);

//functions which are used for the parallel version of the bitonic sort
//imperative openmp version
void impBitonicSortP(void);
//recursive openmp version
void sortP(void);
void bitonicMergeP(int lo, int cnt, int dir);
void recBitonicSortP(int lo, int cnt, int dir);


/** the main program **/ 
int main(int argc, char **argv) {

  if (argc != 3) 
  {
    printf("arguments must be three\n");
    exit(1);
  }
  //check the q and p
  if((atoi(argv[1])<12)||(atoi(argv[1])>24)||(atoi(argv[2])<0)||(atoi(argv[2])>8))
  {
    printf("check again the q and p\n");
    exit(1);
  }
  N = 1<<atoi(argv[1]);//size of the array
  NUM_OF_THREADS=1<< atoi(argv[2]);//2^p threads
  a = (int *) malloc(N * sizeof(int));
  if(a==NULL)
  {
    printf("Not enough memory for the array\n");
    exit(1);
  }

  /* counting the time of qsort algorithm */
  init();
  gettimeofday (&startwtime, NULL);
  qsort(a,N,sizeof(int), cmpfunc_asc);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("QSORT wall clock time = %f\n", seq_time);
  test();

  /* counting the time of the sequential imperative version */
  init();
  gettimeofday (&startwtime, NULL);
  impBitonicSort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Imperative wall clock time = %f\n", seq_time);

  test();

 
  /* counting the time of the imperative OpenMP version */
  init();
  gettimeofday (&startwtime, NULL);
  impBitonicSortP();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Imperative OPENMP clock time = %f\n", seq_time);
  test();

  /* counting the time of the recursive sequential version */
  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive wall clock time = %f\n", seq_time);
  test();

  /* counting the time of the recursive OpenMP version */
  init();
  gettimeofday (&startwtime, NULL);
  sortP();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive OPENMP clock time = %f\n", seq_time);

  test();

  // print();
}

/** -------------- SUB-PROCEDURES  ----------------- **/ 

/** procedure test() : verify sort results **/
void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}


/** procedure init() : initialize array "a" with data **/
void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}


/** INLINE procedure exchange() : pair swap **/
inline void exchange(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}



/** procedure compare() 
   The parameter dir indicates the sorting direction, ASCENDING 
   or DESCENDING; if (a[i] > a[j]) agrees with the direction, 
   then a[i] and a[j] are interchanged.
**/
inline void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j])) 
    exchange(i,j);
}




/** Procedure bitonicMerge() 
   It recursively sorts a bitonic sequence in ascending order, 
   if dir = ASCENDING, and in descending order otherwise. 
   The sequence to be sorted starts at index position lo,
   the parameter cbt is the number of elements to be sorted. 
 **/
void bitonicMerge(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    int i;
    for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
    bitonicMerge(lo, k, dir);
    bitonicMerge(lo+k, k, dir);
  }
}



/** function recBitonicSort() 
    first produces a bitonic sequence by recursively sorting 
    its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order 
 **/
void recBitonicSort(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    recBitonicSort(lo, k, ASCENDING);
    recBitonicSort(lo+k, k, DESCENDING);
    bitonicMerge(lo, cnt, dir);
  }
}


/** function sort() 
   Caller of recBitonicSort for sorting the entire array of length N 
   in ASCENDING order
 **/
void sort() {
  recBitonicSort(0, N, ASCENDING);
}



/*
  imperative version of bitonic sort
*/
void impBitonicSort() {

  int i,j,k;
  
  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=0; i<N; i++) {
	int ij=i^j;
	if ((ij)>i) {
	  if ((i&k)==0 && a[i] > a[ij]) 
	      exchange(i,ij);
	  if ((i&k)!=0 && a[i] < a[ij])
	      exchange(i,ij);
	}
      }
    }
  }
}

// #-----------IMPERATIVE OPENMP VERSION----------------------------------#

//parallel version of the imperative version of the bitonic sort
void impBitonicSortP() {

  int i,j,k;
  
  for (k=2; k<=N; k=2*k) 
  {
    for (j=k>>1; j>0; j=j>>1) 
    {
      //excecuting the for loop in parallel because we can see an independence
      //each parallel excecution has its own i which means that i is private
      //and they have the same a j k as a result a,j,k are shared variables
      #pragma omp parallel for num_threads(NUM_OF_THREADS) private(i) shared(a,j,k)
      for (i=0; i<N; i++) 
      {
	int ij=i^j;
	if ((ij)>i) {
	  if ((i&k)==0 && a[i] > a[ij]) 
	      exchange(i,ij);
	  if ((i&k)!=0 && a[i] < a[ij])
	      exchange(i,ij);
	}
      }
    }
  }
}

// #-----------RECURSIVE OPENMP VERSION----------------------------------#

/*
it is just the same as the bitonicmerge sequential but with different name
because if we make also this function in parallel this has as a result
that the parallel algorithm will be slower than the current parallel
version
*/
void bitonicMergeP(int lo, int cnt, int dir) 
{
  if (cnt>1) 
  {
    int k=cnt/2;
    int i;
    for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
    bitonicMergeP(lo, k, dir);   
    bitonicMergeP(lo+k, k, dir);       
  }
}

/*
We can imagine the recursive version as a tree in which each level has double nodes
in comparison with nodes in previous level.For example we start with one call
of the recursive version,then we can continue with two calls ,the with four calls etc.
So in this version we want to give to all nodes the same number of threads
and thats why we put the if NUM_THREADS/N/cnt > 1 because this command does this 
specific work.
If we cant do further parallelism we continue sequentially but using qsort because it is
faster.
*/

void recBitonicSortP(int lo, int cnt, int dir) 
{
  if (cnt>1) 
  {
    int k=cnt/2;
    //we have to see if we can create threads
    if((NUM_OF_THREADS/(N/cnt))>1)
    {
      #pragma omp parallel sections
      {
        #pragma omp section
        recBitonicSortP(lo, k, ASCENDING);   
        #pragma omp section
        recBitonicSortP(lo+k, k, DESCENDING);
      }
    }
    else//if cant create threads do it sequentially but use qsort instead
    {
      qsort (a+lo,cnt,sizeof(int),dir==ASCENDING?cmpfunc_asc:cmpfunc_dsc);
    }
    bitonicMergeP(lo, cnt, dir);
  }
}

//sortP just calls the parallel recursive version
void sortP() 
{
  recBitonicSortP(0, N, ASCENDING);
}


//sorting in ASCENDING mode
int cmpfunc_asc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}

//if we want to sort in DESCENDING mode
int cmpfunc_dsc (const void * a, const void * b) 
{
   return ( *(int*)b - *(int*)a );
}
