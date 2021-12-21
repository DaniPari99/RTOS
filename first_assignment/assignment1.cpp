/*
--------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------First assignment-----------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------------------------------------
*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/types.h>

// constant value corresponding to the number of tasks
#define NTASKS 4

// number of innerloop and outerloop of the function waste_time()
#define INNERLOOP 1500
#define OUTERLOOP 600

// code of periodic tasks
void task1_code( );
void task2_code( );
void task3_code( );
void task4_code( );

//function for wasting time
void waste_time();


//characteristic function of the thread, only for timing and synchronization
//periodic tasks
void *task1( void *);
void *task2( void *);
void *task3( void *);
void *task4( void *);

//defining some global variable
long int periods[NTASKS];
struct timespec next_arrival_time[NTASKS];
double WCET[NTASKS]; //array of worst case execution time for each task
double max_sched_time[NTASKS];
pthread_attr_t attributes[NTASKS];
pthread_t thread_id[NTASKS];
struct sched_param parameters[NTASKS];
int missed_deadlines[NTASKS]; //array di missed deadlines

//shared variables
int T1T2, T1T4, T2T3;

//blocking time
float B11, B12;
float B21, B22;
float B31;
float B41;



pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER; //initialization of the mutex to protect critical section T1T2
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER; //initialization of the mutex to protect critical section T1T4
pthread_mutex_t mutex3 = PTHREAD_MUTEX_INITIALIZER; //initialization of the mutex to protect critical section T2T3

int main()
{
  	// set task periods in nanoseconds
	//the first task has period 80 millisecond
	//the second task has period 100 millisecond
	//the third task has period 160 millisecond
	//the fourth task has period 200 millisecond
	//you can already order them according to their priority; 
	//if not, you will need to sort them
  	periods[0]= 80000000; //in nanoseconds
  	periods[1]= 100000000; //in nanoseconds
  	periods[2]= 160000000; //in nanoseconds
  	periods[3]= 200000000; //in nanoseconds

	//this is not strictly necessary, but it is convenient to
	//assign a name to the maximum and the minimum priotity in the
	//system. We call them priomin and priomax.
  	struct sched_param priomax;
  	priomax.sched_priority=sched_get_priority_max(SCHED_FIFO);
  	struct sched_param priomin;
  	priomin.sched_priority=sched_get_priority_min(SCHED_FIFO);
  	
  	// set the maximum priority to the current thread (you are required to be
  	// superuser). Check that the main thread is executed with superuser privileges
	// before doing anything else.
 	if (getuid() == 0)	
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax); 

  	// execute all tasks in standalone modality in order to measure execution times
  	// (use gettimeofday).
  	int i;
  	for (i =0; i < NTASKS; i++)
    	{
    		// compute the WCET 40 times for each task and then get the maximum of the 40 values
		WCET[i] = 0;
		for(int j = 0; j < 40; j++)
		{
			// initializa time_1 and time_2 required to read the clock
			struct timespec time_1, time_2;
			clock_gettime(CLOCK_REALTIME, &time_1);

			//we should execute each task more than one for computing the WCET
	 	     	if (i==0)
				task1_code();
      			if (i==1)
				task2_code();
      			if (i==2)
				task3_code();	
			if (i==3)
				task4_code();
      		

			clock_gettime(CLOCK_REALTIME, &time_2); 


			// compute the Worst Case Execution Time
      			max_sched_time[i]= 1000000000*(time_2.tv_sec - time_1.tv_sec) 
			       +(time_2.tv_nsec-time_1.tv_nsec); // tempo già in ns
			       
			//computation of the maximum WCET for each task
			if(max_sched_time[i] > WCET[i]) 
			{
				WCET[i] = max_sched_time[i];
			}
		
			       
		}
      		printf("\nWorst Case Execution Time %d=%f \n", i, WCET[i]);
    	}
    	
    	// compute Ui
    	double U[NTASKS];
    	
    	if(B21 > B41)
    	{
   		U[0] = (WCET[0]/periods[0])+(B21/periods[0]);
 	}
    	else
    	{
    		U[0] = (WCET[0]/periods[0])+(B41/periods[0]);
    	}
    	
    	
    	if(B21 > B41)
    	{
   		U[1] = (WCET[1]/periods[1])+(B21/periods[1]);
 	}
    	else
    	{
    		U[1] = (WCET[1]/periods[1])+(B41/periods[1]);
    	} 
    	
    		
    		U[2] = (WCET[0]/periods[0])+(WCET[1]/periods[1])+(B41/periods[2]);
    		
    		
    		U[3] = (WCET[0]/periods[0])+(WCET[1]/periods[1])+(WCET[2]/periods[2])+(WCET[3]/periods[3]);

	// compute Ulubi
    	double Ulub[NTASKS];
    	int j;
    	for (j =0; j < NTASKS; j++)
    	{ 
    		Ulub[j] = (j+1)*(pow(2.0 , (1.0/(j+1)))-1);
    	}
    		
	//check the sufficient conditions: if they are not satisfied, exit  
  	if (U[0] > Ulub[0] || U[1] > Ulub[1] || U[2] > Ulub[2] || U[3] > Ulub[3])
    	{
      		printf("\n U1=%lf U2=%lf U3=%lf U4=%lf, Ulub1=%lf, Ulub2=%lf, Ulub3=%lf, Ulub4=%lf Non schedulable Task Set\n\n", U[0], U[1], U[2], U[3], Ulub[0], Ulub[1], Ulub[2], Ulub[3]);
      		return(-1);
    	}
  	printf("\n U1=%lf U2=%lf U3=%lf U4=%lf, Ulub1=%lf, Ulub2=%lf, Ulub3=%lf, Ulub4=%lf Schedulable Task Set\n\n", U[0], U[1], U[2], U[3], Ulub[0], Ulub[1], Ulub[2], Ulub[3]);
  	fflush(stdout);
  	sleep(5);

       // set the minimum priority to the current thread: this is now required because 
	//we will assign higher priorities to periodic threads to be soon created
	//pthread_setschedparam
  	if (getuid() == 0)	//ora assegno la priorità minore al main thread
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);
    		
    	// set the attributes of each task, including scheduling policy and priority
  	for (i =0; i < NTASKS; i++)
    	{
		//initializa the attribute structure of task i
      		pthread_attr_init(&(attributes[i]));

		//set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		//not inherited from the main thread (pthread_attr_setinheritsched) 
      		pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
		pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO);

		//properly set the parameters to assign the priority inversely proportional 
		//to the period
      		parameters[i].sched_priority = priomin.sched_priority+NTASKS - i;

		//set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
      		pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    	}
    	
    	//set the attributes of each mutex
    	pthread_mutexattr_t mymutexattr;
    	pthread_mutexattr_init(&mymutexattr);
    	pthread_mutexattr_setprotocol(&mymutexattr, PTHREAD_PRIO_PROTECT);
    	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[0].sched_priority);
    	pthread_mutex_init(&mutex1, &mymutexattr);
    	pthread_mutex_init(&mutex2, &mymutexattr);
    	pthread_mutexattr_setprioceiling(&mymutexattr, parameters[1].sched_priority);
    	pthread_mutex_init(&mutex3, &mymutexattr);
    	
    	
    	//declare the variable to contain the return values of pthread_create	
  	int iret[NTASKS];

	//declare variables to read the current time
	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1);
	
	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (i = 0; i < NTASKS; i++)
    	{
		long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i];
		
		//then we compute the end of the first period and beginning of the next one
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000;
       		missed_deadlines[i] = 0;
    	}
    	
    	// create all threads
  	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
  	iret[3] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);

  	// join all threads
  	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
  	pthread_join( thread_id[3], NULL);
  	
  	// print the missed deadlines for each task
  	for (i = 0; i < NTASKS; i++)
    	{
      		printf ("\nMissed Deadlines Task %d=%d\n\n", i, missed_deadlines[i]);
		fflush(stdout);
    	}
  	exit(0);
}

// task1: shall write something into T1T2 and into T1T4 	
void task1_code()
{
	printf(" 1[ "); fflush(stdout);
  	struct timespec time1, time2, time3, time4; // objects time to calculate the blocking time Bij
  	int a = 2;
  	int b = 3;
  		
  	pthread_mutex_lock(&mutex1);	// function to get the mutex1
  	clock_gettime(CLOCK_REALTIME, &time1);
  	T1T2 = a*b;	//writing in T1T2
  	
  	//waste_time();
  	//waste_time();
  	//waste_time();
  	
  	clock_gettime(CLOCK_REALTIME, &time2);
  	pthread_mutex_unlock(&mutex1);	// function to release the mutex1
  	B11 = 1000000000*(time2.tv_sec-time1.tv_sec) + (time2.tv_nsec-time1.tv_nsec);	// computation of Bij
  		
  	pthread_mutex_lock(&mutex2);
  	clock_gettime(CLOCK_REALTIME, &time3);
  	T1T4 = a+b;
  	clock_gettime(CLOCK_REALTIME, &time4);
  	pthread_mutex_unlock(&mutex2);
  	B12 = 1000000000*(time2.tv_sec-time1.tv_sec) + (time2.tv_nsec-time1.tv_nsec); 
  	printf(" ]1 "); fflush(stdout); 		
}
  	
  	//thread code for task_1 (used only for temporization)
void *task1( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

    	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		// execute application specific code
		task1_code();
		
		//check if there are missed deadlines and in the case that there are increment the counter
		struct timespec missed_time;
		clock_gettime(CLOCK_REALTIME, &missed_time);
		if((1000000000*(missed_time.tv_sec) + missed_time.tv_nsec) > (1000000000*(next_arrival_time[0].tv_sec) + next_arrival_time[0].tv_nsec))
		{
			missed_deadlines[0] += 1;
		}
		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL);

		long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
		next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
		

		
    	}
}
 
// task2: shall read something from T1T2 and write something into T2T3
void task2_code()
{
	printf(" 2[ "); fflush(stdout);
  	struct timespec time1, time2, time3, time4; // objects time to calculate the blocking time Bij
  	int a;
  	int b = 1;

  	pthread_mutex_lock(&mutex1);	// function to get the mutex1
  	clock_gettime(CLOCK_REALTIME, &time1);
  	a = T1T2*2;	//reading from T1T2
  	
  	//waste_time();
  	//waste_time();
  	
  	clock_gettime(CLOCK_REALTIME, &time2);
  	pthread_mutex_unlock(&mutex1);  	// function to release the mutex1
  	//pthread_mutex_unlock(&mutex3);		
  	B21 = 1000000000*(time2.tv_sec-time1.tv_sec) + (time2.tv_nsec-time1.tv_nsec);	// computation of Bij
  		  			  		
  	pthread_mutex_lock(&mutex3);	// function to get the mutex3
  	clock_gettime(CLOCK_REALTIME, &time3);
  	T2T3 = a+b;	//writing into T2T3
  	clock_gettime(CLOCK_REALTIME, &time4);
  	pthread_mutex_unlock(&mutex3);  	// function to release the mutex3	
  	B22 = 1000000000*(time2.tv_sec-time1.tv_sec) + (time2.tv_nsec-time1.tv_nsec);	// computation of Bij
  	printf(" ]2 "); fflush(stdout);
  		
}
  	
void *task2( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

    	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		// execute application specific code
		task2_code();
		
		//check if there are missed deadlines and in the case that there are increment the counter
		struct timespec missed_time;
		clock_gettime(CLOCK_REALTIME, &missed_time);
		if((1000000000*(missed_time.tv_sec) + missed_time.tv_nsec) > (1000000000*(next_arrival_time[1].tv_sec) + next_arrival_time[1].tv_nsec))
		{
			missed_deadlines[1] += 1;
		}

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);

		long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
		next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;
	
    	}
}

// task3: shall read something from T2T3
 void task3_code()
{
	printf(" 3[ "); fflush(stdout);
  	struct timespec time1, time2; 	// objects time to calculate the blocking time Bij
  	int a;
  		
  	pthread_mutex_lock(&mutex3);		// function to get mutex3		
  	clock_gettime(CLOCK_REALTIME, &time1);
  	a = T2T3*2;
  	
  	//waste_time();
  	//waste_time();
  	
  	clock_gettime(CLOCK_REALTIME, &time2);
  	pthread_mutex_unlock(&mutex3);  	// function to release mutex3	
  	B31 = 1000000000*(time2.tv_sec-time1.tv_sec) + (time2.tv_nsec-time1.tv_nsec);		// computation of Bij
  	printf(" ]3 "); fflush(stdout);

}

void *task3( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

    	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		// execute application specific code
		task3_code();
		
		//check if there are missed deadlines and in the case that there are increment the counter
		struct timespec missed_time;
		clock_gettime(CLOCK_REALTIME, &missed_time);
		if((1000000000*(missed_time.tv_sec) + missed_time.tv_nsec) > (1000000000*(next_arrival_time[2].tv_sec) + next_arrival_time[2].tv_nsec))
		{
			missed_deadlines[2] += 1;
		}

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);

		long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
		next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
		
		
    	}
}


// task4: shall read something from T1T4
 void task4_code()
{
	printf(" 4[ "); fflush(stdout);
  	struct timespec time1, time2; 	// objects time to calculate the blocking time Bij
  	int a;
  		
  	pthread_mutex_lock(&mutex2);		// function to get mutex2
  	clock_gettime(CLOCK_REALTIME, &time1);
  	a = T1T4*3;
  	
  	//waste_time();
  	//waste_time();
  	
  	clock_gettime(CLOCK_REALTIME, &time2);
  	pthread_mutex_unlock(&mutex2);	// function to release mutex2  		
  	B41 = 1000000000*(time2.tv_sec-time1.tv_sec) + (time2.tv_nsec-time1.tv_nsec);	// computation of Bij
  	printf(" ]4 "); fflush(stdout);
}

void *task4( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

    	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		// execute application specific code
		task4_code();
		
		//check if there are missed deadlines and in the case that there are increment the counter
		struct timespec missed_time;
		clock_gettime(CLOCK_REALTIME, &missed_time);
		if((1000000000*(missed_time.tv_sec) + missed_time.tv_nsec) > (1000000000*(next_arrival_time[3].tv_sec) + next_arrival_time[3].tv_nsec))
		{
			missed_deadlines[3] += 1;
		}

		clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);

		long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
		next_arrival_time[3].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[3].tv_sec= next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
		
		
    	}
}

// function to waste time
void waste_time()
{
	for (int i = 0; i < OUTERLOOP; i++)
    	{
      		for (int j = 0; j < INNERLOOP; j++);		
			double uno = rand()*rand()%10;
    	}
}
