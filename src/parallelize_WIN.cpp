/*  P A R A L L E L I Z E

    Some basic support to thread parallelization.

    Davide Bucci, 2014-2016
    Lionel Bastard, 2026

*/

/* This file is part of cRCWA.

    cRCWA is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.
*/

#include "parallelize.h"
#include "compat.h"
#include <iostream>

static bool windows_parallel_warning_shown = false;

void parallelize(struct thread_data *td, void* (pfunc)(void*), int instance,
    int nmax, int number_of_instances)
{
    if (!windows_parallel_warning_shown && nmax>1) {
        std::cerr << "[parallelize] Warning: parallel execution not available "
            "on Windows build; running serially." << std::endl;
        windows_parallel_warning_shown = true;
    }
    td->close_process = false;
    pfunc((void*)td);
}

void *parallelize_finish(struct thread_data *my_data)
{
    return NULL;
}

void waitSemaphoreIO(void) { /* no-op */ }
void postSemaphoreIO(void) { /* no-op */ }

void wait_for_threads(vector<struct thread_data> td, int number_of_instances)
{ /* no-op */ }
