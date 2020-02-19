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
#include "cilk/cilk.h"
#include "cilk/cilk_api.h"
#include <math.h>

struct timeval startwtime, endwtime;
double seq_time;


int N;          // data array size
int *a;         // data array to be sorted

const int ASCENDING  = 1;
const int DESCENDING = 0;

int NUM_OF_THREADS;  // number of threads

void init(void);
void print(void);
void sort(void);
void test(void);
inline void exchange(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSort(int lo, int cnt, int dir);
void impBitonicSort(void);

//function useful for the use of quicksort
int cmpfunc (const void * a, const void * b);

// functions which are made for the parallel version
void sortP();
void impBitonicSortP();
void bitonicMergeP(int lo, int cnt, int dir);
void recBitonicSortP(int lo, int cnt, int dir);

/** the main program **/ 
int main(int argc, char **argv) 
{

  if (argc != 3) // if args are not 3 then do not continue
  {
    printf("check the number of arguments (must be 3)\n");
    exit(1);
  }
  //argv[1]=q  argv[2]=p
  //checking the q and p if they are between the limits
  if((atoi(argv[1])<12)||(atoi(argv[1])>24)||(atoi(argv[2])<0)||(atoi(argv[2])>8))
  {
    printf("check again the q and p \n");
    exit(1);
  }
  N = 1<<atoi(argv[1]);  // N =2^q elements
  a = (int *) malloc(N * sizeof(int));
  //check if there is enough space for the array
  if(a==NULL)
  {
    printf("not enough memory for the array\n");
    exit(1);
  }
  NUM_OF_THREADS=1<<atoi(argv[2]); // num of threads = 2^p
  /* counting the time of qsort algorithm */
  init();
  gettimeofday (&startwtime, NULL);
  qsort(a,N,sizeof(int), cmpfunc);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("QSORT wall clock time = %f\n", seq_time);
  test();
  /* counting the time of the recursive sequential version */
  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive Sequential wall clock time = %f\n", seq_time);
  test();

  /* counting the time of the recursive parallel version */
  init();
  gettimeofday (&startwtime, NULL);
  sortP();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive PARARELL(cilk) clock time = %f\n", seq_time);
  test();
  /* counting the time of the sequential imperative version */
  init();
  gettimeofday (&startwtime, NULL);
  impBitonicSort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("imperative Sequential wall clock time = %f\n", seq_time);
  test();

  /* counting the time of the imperative CilkPlus version */
  init();
  gettimeofday (&startwtime, NULL);
  impBitonicSortP();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("imperative PARARELL(cilk) wall clock time = %f\n", seq_time);

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


//------------------------CilkPlus Versions--------------------------------


void sortP() 
{
  //convert the Num_of_threads from integer into string
  char THREADS[21];
  sprintf(THREADS, "%d",NUM_OF_THREADS);
  __cilkrts_set_param ("nWorkers",THREADS);
  recBitonicSortP(0, N, ASCENDING);
}




void impBitonicSortP() 
{
  //convert the Num_of_threads from ineteger into string
  char THREADS[21];
  sprintf(THREADS, "%d",NUM_OF_THREADS);
  __cilkrts_set_param ("nWorkers",THREADS);
  int i,j,k; 
  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) 
    {
      //we are excecuting the for loop in parallel
      //because we can see the independence
      cilk_for (i=0; i<N; i++) 
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

/*
  we do not touch the bitonicMerge with Cilk
  because this will make the parallelism slower
  than it is now.
  this is a result from running this program many times
  and counting times
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


void recBitonicSortP(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    
    cilk_spawn recBitonicSortP(lo, k, ASCENDING);
    recBitonicSortP(lo+k, k, DESCENDING);
    cilk_sync;//wait here for the recbitonic sort to finish
    //call always the merge
    bitonicMergeP(lo, cnt, dir);
  }
}

//sorting in ASCENDING mode with this function(qsort)
int cmpfunc (const void * a, const void * b)
{
   return ( *(int*)a - *(int*)b );
}






