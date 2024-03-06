#ifndef RINGBUFFER_H
#define RINGBUFFER_H

#include "main.h"
#include <vector>

using namespace std;

class RingBuffer
{
public:
    ~RingBuffer();

    void initialize(int size);

    void write(const int16_t data);
    void write(const int16_t *data, const int count);
    
    int read(int16_t *data, const int count);

    bool isFull();
    bool isDataAvailable();
    int bufferSize();
    int maximumSize();
    void printData();

private:
    vector<int16_t> buffer;
    int size, head, tail;
    bool isEmpty;

    int16_t read();
};

#endif