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
#include <pthread.h>

struct timeval startwtime, endwtime;
double seq_time;


int N;          // data array size
int *a;         // data array to be sorted

const int ASCENDING  = 1;
const int DESCENDING = 0;

int NUM_OF_THREADS;//2^p threads (maximum)
int activethreads;//threads that are active at a specific time

pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;//useful for the change of the variable: activethreads



void init(void);
void print(void);
void sortP(void);//function for the parallel algorithm
void test(void);
inline void exchange(int i, int j);
void compare(int i, int j, int dir);

/* functions for the parallelism */
void bitonicMergeP(int lo, int cnt, int dir);
void recBitonicSortP(int lo, int cnt, int dir);

/* functions for the call of qsort */
int cmpfunc_asc (const void * a, const void * b);//for ascending
int cmpfunc_dsc (const void * a, const void * b);//for descending

/* next 3 functions for the sequential algorithm */
void recBitonicSort(int lo, int cnt, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void sort(void);

/* functions called by created threads */
void *recBitonicThread (void* arg);
void *MergeBitonicThread (void* arg);

typedef struct
{ 
  int lo,cnt,dir;
}thread_struct;


/** the main program **/ 
int main(int argc, char **argv) 
{
  //check if args are 3
  if (argc != 3) 
  {
    printf("check the number of arguments (must be 3)\n");
    exit(1);
  }
  //check the q and p
  if((atoi(argv[1])<12)||(atoi(argv[1])>24)||(atoi(argv[2])<0)||(atoi(argv[2])>8))
  {
    printf("check again the q and p\n");
    exit(1);
  }

  N = 1<<atoi(argv[1]);//size of the array
  a = (int *) malloc(N * sizeof(int));
  if(a==NULL)
  {
    printf("not enough memory for the array\n");
    exit(1);
  }
  NUM_OF_THREADS=1<<atoi(argv[2]);//2^p threads
  /* counting the time of qsort algorithm */
  init();
  gettimeofday (&startwtime, NULL);
  qsort(a,N,sizeof(int),cmpfunc_asc);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Qsort wall clock time = %f\n", seq_time);
  test();
  /* counting the time of the recursive sequential version */
  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive SERIAL clock time = %f\n", seq_time);
  test();
  /* counting the time of the recursive parallel version */
  init();
  gettimeofday (&startwtime, NULL);
  sortP();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive PARALLEL(THREADS) clock time = %f\n", seq_time);
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


void recBitonicSortP(int lo, int cnt, int dir) 
{
  if (cnt>1) 
  {
    int k=cnt/2;
    //if we can create more threads do the algorithm in parallel
    if(activethreads<NUM_OF_THREADS)
    {
      //thread1 will "call" the recBitonicSortP with dir=ASCENDING
      //thread1 will "call" the recBitonicSortP with dir=DESCENDING
      pthread_t thread1,thread2;
      // struct1 for ASCENDING struct2 for DESCENDING
      thread_struct struct1,struct2;
      //we use mutex because we want only one thread at a time to write in this variable
      //we dont want race conditions
      pthread_mutex_lock (&mutex1);
      activethreads+=2;
      pthread_mutex_unlock(&mutex1);

      struct1.lo=lo;
      struct1.cnt=k;
      struct1.dir=ASCENDING;

      struct2.lo=lo+k;
      struct2.cnt=k;
      struct2.dir=DESCENDING;
      //creating the two threads
      pthread_create (&thread1,NULL,recBitonicThread,&struct1);
      pthread_create (&thread2,NULL,recBitonicThread,&struct2);
      //waiting for threads to end their job
      pthread_join ( thread1, NULL );
      pthread_join ( thread2, NULL );
      //now we can use mutex again because we know that 2 threads are "dead"
      pthread_mutex_lock (&mutex1);
      activethreads-=2;
      pthread_mutex_unlock(&mutex1);
    }
    //if we cant create more threads do the algorithm sequentially
    //but use qsort instead of the recursive version
    //because qsort is quickier
    else
    {
      qsort (a+lo,cnt,sizeof(int),dir==ASCENDING?cmpfunc_asc:cmpfunc_dsc);
    }
    //calling the merge function always
    bitonicMergeP(lo, cnt, dir);
  }
}



void bitonicMergeP(int lo, int cnt, int dir) 
{
  if (cnt>1) 
  {
    int k=cnt/2;
    int i;
    //we are doing the same job as at the function recBitonicSortP
    //see the comments there
    if(activethreads<NUM_OF_THREADS)
    {
      
      pthread_t thread1,thread2;
      thread_struct struct1,struct2;

      pthread_mutex_lock (&mutex1);
      activethreads+=2;
      pthread_mutex_unlock(&mutex1);

      struct1.lo=lo;
      struct1.cnt=k;
      struct1.dir=dir;

      struct2.lo=lo+k;
      struct2.cnt=k;
      struct2.dir=dir;
      
      for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
      
      pthread_create (&thread1,NULL,MergeBitonicThread,&struct1);
      pthread_create (&thread2,NULL,MergeBitonicThread,&struct2);
      pthread_join ( thread1, NULL );
      pthread_join ( thread2, NULL );
      pthread_mutex_lock (&mutex1);
      activethreads-=2;
      pthread_mutex_unlock(&mutex1);
          
    }
    else//do the algorithm sequentially
    {
      for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
      bitonicMergeP(lo, k, dir);
      bitonicMergeP(lo+k, k, dir);
    }
    
  }
}


/** function sort() 
   Caller of recBitonicSort for sorting the entire array of length N 
   in ASCENDING order
 **/
void sortP() 
{
  activethreads = 0;
  recBitonicSortP(0, N, ASCENDING);
}

//if we want to sort in ASCENDING mode
int cmpfunc_asc (const void * a, const void * b) 
{
   return ( *(int*)a - *(int*)b );
}
//if we want to sort in DESCENDING mode
int cmpfunc_dsc (const void * a, const void * b) 
{
   return ( *(int*)b - *(int*)a );
}

//here we are just taking the lo,cnt and dir from the arg
//and calling the recBitonicSortP with these variables as arguments
void *recBitonicThread (void* arg)
{
  if ( ((thread_struct *)arg)-> cnt > 1 ) 
  {
    recBitonicSortP( ((thread_struct *)arg)-> lo, ((thread_struct *)arg)-> cnt, ((thread_struct *)arg)-> dir );
  }
}

//the same as void *recBitonicThread (void* arg)
//but we calling the bitonicMergeP this time
void *MergeBitonicThread (void* arg)
{
  if ( ((thread_struct *)arg)-> cnt > 1 ) 
  {
    bitonicMergeP( ((thread_struct *)arg)-> lo, ((thread_struct *)arg)-> cnt, ((thread_struct *)arg)-> dir );
  }
}



/* ----------------------------Sequential functions---------------------------------------------- */

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

void sort() 
{
  recBitonicSort(0, N, ASCENDING);
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

