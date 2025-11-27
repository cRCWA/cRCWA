#ifndef PARALLELIZE_H
#define PARALLELIZE_H

#include "structure.h"

struct thread_data{
   int  instance_number;
   structure *p;
   void *payload;
   bool close_process;
   int *thread_counter;
   pthread_t thread;
};

void parallelize(struct thread_data *td, void* (pfunc)(void*), int i, int nmax,
	int number_of_instances);

void *parallelize_finish(struct thread_data *td);
void wait_for_threads(struct thread_data *td, int number_of_instances);

void waitSemaphoreIO(void);
void postSemaphoreIO(void);

#endif