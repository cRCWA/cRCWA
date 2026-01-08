/*  P A R A L L E L I Z E

    Some basic support to thread parallelization.

    Davide Bucci, 2014-2016

*/

/* This file is part of cRCWA.

    cRCWA is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.
    
    cRCWA is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
    details.
    
    You should have received a copy of the GNU General Public License along
    with cRCWA. If not, see <https://www.gnu.org/licenses/>. 

    Davide Bucci, 2008-2026
    Jérôme Michallon, 2012-2014
*/

#include <pthread.h>
#include <unistd.h>

#include <semaphore.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <stdlib.h>

#include "parallelize.h"

sem_t *mutex_thread_counter; // Semaphore for thread counter access.
#define BUFFER_SIZE 255
char semaphore_thread_name[BUFFER_SIZE+1];

sem_t *mutex_io_access; // Semaphore for I/O.
char semaphore_io_access_name[BUFFER_SIZE+1];


/** Main scheduler function: launch a new thread with the given function
    with the given data, up to a certain number of threads. This function
    returns immediately when the number of created threads is less than the
    maximum number of allowable threads.
    Synchronisation is to be done with wait_for_threads().
*/
void parallelize(struct thread_data *td, void* (pfunc)(void*), int instance,
    unsigned int nmax, unsigned int number_of_instances)
{
    int rc;
    pthread_t thread;
    int p=getpid();

    if(nmax==0) {
        cerr<<"Warning: parallelization scheduler can not be set with a"
            " maximum number of instances equal to zero."<<endl;
    }

    if(mutex_thread_counter==NULL && nmax>1) {
        snprintf(semaphore_thread_name,BUFFER_SIZE,
            "/semaphore_thread_c_%d", p);
        mutex_thread_counter=sem_open(semaphore_thread_name, O_CREAT, 0644, 1);
        if(mutex_thread_counter==SEM_FAILED) {
            cerr<<"Could not set up a semaphore (thread counter), and "
            "I really need it!";
            cerr<<endl<<"errno="<<errno<<endl;
            cout<<"Thread counter semaphore not set."<<endl;
            exit(-1);
        }
    }
    if(mutex_io_access==NULL && nmax>1) {
        snprintf(semaphore_io_access_name,BUFFER_SIZE,
            "/semaphore_io_access_%d", p);
        mutex_io_access=sem_open(semaphore_io_access_name, O_CREAT, 0644, 1);
        if(mutex_io_access==SEM_FAILED) {
            cerr<<"Could not set up a semaphore (IO access), and "
            "I really need it!";
            cerr<<endl<<"errno="<<errno<<endl;
            cout<<"Thread counter semaphore not set."<<endl;
            exit(-1);
        }
    }

    if(*td->thread_counter >= nmax || instance==number_of_instances-1) {
        td->close_process=false;
        td->thread=pthread_self();
        pfunc((void *)td);
    } else {
        sem_wait (mutex_thread_counter);
        ++(*td->thread_counter);
        sem_post (mutex_thread_counter);

        td->close_process=true;
        rc = pthread_create(&thread, NULL, pfunc, (void *)td);
        waitSemaphoreIO();
        //printf("-----> Thread created.\n");
        postSemaphoreIO();
        td->thread=thread;
    }


}

/** Handle the end of a parallelised function.

*/
void *parallelize_finish(struct thread_data *my_data)
{

    if(my_data->close_process) {
        sem_wait (mutex_thread_counter);
        --*(my_data->thread_counter);
        sem_post (mutex_thread_counter);
        pthread_exit(NULL);
    } else {
        return NULL;
    }
    return NULL;
}

void waitSemaphoreIO(void)
{
    if(mutex_io_access!=NULL) sem_wait (mutex_io_access);
}


void postSemaphoreIO(void)
{
    if(mutex_io_access!=NULL) sem_post (mutex_io_access);
}

/**
    Wait until all launched threads have finished their tasks.
 */
void wait_for_threads(struct thread_data *td, int number_of_instances)
{
    bool to_cont=true;
    waitSemaphoreIO();
    //printf("Synchronization\n");
    postSemaphoreIO();
    void *retval;

    for(int i=0; i<number_of_instances; ++i) {
         if(pthread_equal(pthread_self(),td[i].thread)==0 &&
            pthread_join(td[i].thread, &retval)!=0) {
            waitSemaphoreIO();
            printf("Error! can not join threads.\n");
            postSemaphoreIO();
            break;
         } else {
            waitSemaphoreIO();
            //printf("Thread joined\n");
            postSemaphoreIO();
         }
    }

    // close semaphores
    if(mutex_thread_counter!=NULL) {
        sem_close(mutex_thread_counter);
        sem_unlink(semaphore_thread_name);
        mutex_thread_counter=NULL;
    }

    if(mutex_io_access!=NULL) {
        sem_close(mutex_io_access);
        sem_unlink(semaphore_io_access_name);
        mutex_io_access=NULL;
    }
}
