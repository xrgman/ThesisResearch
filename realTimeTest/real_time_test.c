#include <stdio.h>
#include <time.h>
#include <unistd.h>

#define INTERVAL_MS 100 // Time interval between messages in milliseconds

void print_current_time() {
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    printf("Current time: %ld.%09ld\n", now.tv_sec, now.tv_nsec);
}

int main() {
    struct timespec interval;
    interval.tv_sec = INTERVAL_MS / 1000;
    interval.tv_nsec = (INTERVAL_MS % 1000) * 1000000;

    printf("Starting real-time task...\n");

    while (1) {
        print_current_time();
        nanosleep(&interval, NULL); // Sleep for the specified interval
    }

    return 0;
}