/* Minimal portability helpers for Windows vs POSIX */
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

static inline sem_t sem_open(const char* /*name*/, int /*oflag*/, int /*mode*/, unsigned int /*value*/) {
    return new sem_t_wrapper();
}
static inline int sem_close(sem_t s) { if (s) delete s; return 0; }
static inline int sem_unlink(const char* /*name*/) { return 0; }
static inline int sem_wait(sem_t s) { if (s) { s->m.lock(); return 0; } return -1; }
static inline int sem_post(sem_t s) { if (s) { s->m.unlock(); return 0; } return -1; }
static inline int sem_init(sem_t* s, int /*pshared*/, unsigned int /*value*/) { *s = new sem_t_wrapper(); return 0; }
static inline int sem_destroy(sem_t* s) { if (*s) delete *s; *s = nullptr; return 0; }

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

/* simple getpass replacement: reads password without echoing */
static inline char *getpass(const char *prompt)
{
    static char buf[256];
    if (prompt) fprintf(stderr, "%s", prompt);
    int i = 0;
    int c;
    while ((c = _getch()) != '\r' && c != EOF && i < (int)sizeof(buf) - 1) {
        if (c == '\b') {
            if (i > 0) { --i; }
        } else {
            buf[i++] = (char)c;
        }
    }
    buf[i] = '\0';
    if (prompt) fprintf(stderr, "\n");
    return buf;
}

#else /* POSIX */
#include <unistd.h>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#endif

#endif /* COMPAT_H */
