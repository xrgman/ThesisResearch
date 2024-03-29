#ifndef RINGBUFFER_H
#define RINGBUFFER_H

#include "main.h"
#include <vector>
#include <chrono>

using namespace std;

class RingBuffer
{
public:
    ~RingBuffer();

    void initialize(int size);
    void reset();

    void write(const int16_t data);
    void write(const int16_t data, chrono::time_point<chrono::high_resolution_clock> receivedTime);
    void write(const int16_t *data, const int count);

    int16_t read();
    int16_t read(std::chrono::time_point<std::chrono::high_resolution_clock> &receivedTime);
    int read(int16_t *data, const int count);

    bool isFull();
    bool isDataAvailable();
    chrono::time_point<chrono::high_resolution_clock> getEmptyTime();
    int bufferSize();
    int maximumSize();
    void printStats();
    void printData();

private:
    vector<int16_t> buffer;
    vector<pair<int, chrono::time_point<chrono::high_resolution_clock>>> receivedTimes;
    int size, head, tail;
    bool isEmpty;
    chrono::time_point<chrono::high_resolution_clock> isEmptyTime;

    bool receivedTimeSeen(const std::chrono::time_point<std::chrono::high_resolution_clock> &receivedTime);
};

#endif