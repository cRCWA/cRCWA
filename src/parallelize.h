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

    Davide Bucci, 2008-2025
    Jérôme Michallon, 2012-2014
*/

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