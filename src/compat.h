/* Minimal portability helpers for Windows vs POSIX */


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

    Lionel Bastard, 2025-2026
    Davide Bucci, 2026
*/


#ifndef COMPAT_H
#define COMPAT_H

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <process.h>
#include <conio.h>
#include <cstdio>
#include <mutex>
#include <map>
#include <string>

/* Minimal semaphore shim for Windows builds: maps named/unnamed semaphores
   to a simple mutex. This is sufficient for serializing FFTW/IO calls in
   single-process builds. */
struct sem_t_wrapper { std::mutex m; };
typedef sem_t_wrapper* sem_t;

static inline sem_t sem_open(const char* /*name*/, int /*oflag*/, int /*mode*/,
    unsigned int /*value*/)
{
    return new sem_t_wrapper();
}
static inline int sem_close(sem_t s) { if (s) delete s; return 0; }
static inline int sem_unlink(const char* /*name*/) { return 0; }
static inline int sem_wait(sem_t s)
{
    if (s) { s->m.lock(); return 0; } return -1; 
}
static inline int sem_post(sem_t s)
{
    if (s) { s->m.unlock(); return 0; } return -1;
}
static inline int sem_init(sem_t* s, int /*pshared*/, unsigned int /*value*/) 
{
    *s = new sem_t_wrapper(); return 0;
}

static inline int sem_destroy(sem_t* s)
{
    if (*s) delete *s; *s = nullptr; return 0; 
}

#define SEM_FAILED ((sem_t)NULL)
#ifndef O_CREAT
#define O_CREAT 0
#endif

static inline int getpid_compat() { return _getpid(); }
#define getpid getpid_compat

static inline void usleep(unsigned int usec) { Sleep((usec) / 1000); }

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_E
#define M_E 2.71828182845904523536
#endif

#else /* POSIX */
#include <unistd.h>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

#endif /* COMPAT_H */
