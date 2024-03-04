#ifndef RINGBUFFER_H
#define RINGBUFFER_H

#include "main.h"
#include <vector>

using namespace std;

class RingBuffer
{
public:
    ~RingBuffer();

    void Initialize(int size);

    void Write(const int16_t data);
    void Write(const int16_t *data, const int count);
    int16_t Read();
    int Read(int16_t *data, const int count);

    bool isDataAvailable();
    int bufferSize();
    void printData();

private:
    vector<int16_t> buffer;
    int size, head, tail;
    bool isEmpty;
};

#endif