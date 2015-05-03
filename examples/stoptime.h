#include <time.h>
#include <sys/time.h>


double stoptime(void)
{
    struct timeval	t;
    gettimeofday(&t,NULL);
    return (double) t.tv_sec + t.tv_usec/1000000.0;
}


void exactsleep (double length)
{
    struct timespec requested_time, remaining;
    requested_time.tv_sec=(int) length;
    requested_time.tv_nsec=(int) ((length-(int)length)*1e9);
    nanosleep (&requested_time, &remaining);
}


#ifdef __MMX__
//#include <unistd.h>
#include <stdint.h>
// using the Time Stamp Counter register (64-bit) present since Pentium
// tsc is incremented every cpu tick
// instruction: rdtsc: reads tsc into edx:eax
extern "C" {
    __inline__ uint64_t rdtsc(void)
    {
        uint32_t lo, hi;
        __asm__ __volatile__ (      // serialize
            "xorl %%eax,%%eax \n        cpuid"
            ::: "%rax", "%rbx", "%rcx", "%rdx");
        /* We cannot use "=A", since this would use %rax on x86_64 and return only the lower 32bits of the TSC */
        __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
        return (uint64_t)hi << 32 | lo;
    }
}
#endif
