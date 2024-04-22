#include "ringBuffer.h"
#include "util.h"
#include <iostream>

/// @brief Deconstructor.
RingBuffer::~RingBuffer()
{
    // delete[] buffer;
}

/// @brief Initialize the buffer.
/// @param size Size of the buffer.
void RingBuffer::initialize(int size)
{
    this->size = size;
    this->head = 0;
    this->tail = 0;
    this->isEmpty = true;

    // Allocating the buffer with the right size:
    // buffer = new int16_t[size];
    buffer.resize(size);
    buffer.clear();

    receivedTimes.clear();
    receivedTimes.reserve(20);

    fillArrayWithZeros(buffer.data(), size);
}

/// @brief Reset the buffer.
void RingBuffer::reset()
{
    this->head = 0;
    this->tail = 0;
    this->isEmpty = true;

    buffer.clear();
    receivedTimes.clear();
}

/// @brief Write a single element to the buffer.
/// @param data Data element to write to the buffer.
void RingBuffer::write(const int16_t data)
{
    // Checking for overflow:
    if (head == tail && !isEmpty)
    {
        //spdlog::error("Ringbuffer: Buffer overflow, overwriting oldest data.");

        tail = (tail + 1) % size;
    }

    isEmpty = false;

    // Writing data to front of buffer:
    buffer[head] = data;

    // Increasing write index:
    head = (head + 1) % size;
}

/// @brief Writing an entire data array to the buffer.
/// @param data Array containing the data to be written.
/// @param count Number of elements in the array.
void RingBuffer::write(const int16_t *data, const int count)
{
    for (int i = 0; i < count; ++i)
    {
        write(data[i]);
    }
}

/// @brief Write a single element to the buffer.
/// @param data Data element to write to the buffer.
/// @param receivedTime The time that the batch containing this data was received.
void RingBuffer::write(const int16_t data, chrono::time_point<chrono::high_resolution_clock> receivedTime)
{
    // Saving the received time at current header position:
    if (!receivedTimeSeen(receivedTime))
    {
        receivedTimes.push_back(pair<int, chrono::time_point<chrono::high_resolution_clock>>(head, receivedTime));
    }

    // Performing actual writing to buffer:
    write(data);
}

/// @brief Read a single element from the buffer.
/// @return The read element.
int16_t RingBuffer::read()
{
    // Checking if buffer is not empty:
    if (isEmpty)
    {
        return -1;
    }

    int16_t element = buffer[tail];

    // Increasing read index:
    tail = (tail + 1) % size;

    // If buffer is empty, declare it that way :)
    if (tail == head)
    {
        isEmpty = true;

        // Log time that buffer went empty:
        isEmptyTime = chrono::high_resolution_clock::now();
    }

    return element;
}

int16_t RingBuffer::read(std::chrono::time_point<std::chrono::high_resolution_clock> &receivedTime)
{
    // Grabbing time:
    if (receivedTimes.size() > 0)
    {
        while (receivedTime.time_since_epoch().count() == 0)
        {
            pair<int, chrono::time_point<chrono::high_resolution_clock>> receivedTimePair = receivedTimes[0];

            // Checking if received time can be removed:
            if (receivedTimes.size() > 1)
            {
                pair<int, chrono::time_point<chrono::high_resolution_clock>> receivedTimePair2 = receivedTimes[1];

                if (receivedTimePair2.first < tail)
                {
                    receivedTimes.erase(receivedTimes.begin());

                    continue;
                }
            }

            receivedTime = receivedTimePair.second;
        }
    }

    return read();
}

/// @brief Read multiple elements from the buffer.
/// @param data Array to store the read elements in.
/// @param count Number of elements to read.
/// @return Number of elements that are actually read.
int RingBuffer::read(int16_t *data, const int count)
{
    int readCnt = 0;

    for (int i = 0; i < count; ++i)
    {
        // Checking if buffer is empty:
        if (isEmpty)
        {
            break;
        }

        data[i] = read();

        // Keeping track of the total amount of read elements:
        readCnt++;
    }

    return readCnt;
}

/// @brief Function to check whether the buffer is full.
/// @return Whether or not the buffer is at capacity.
bool RingBuffer::isFull()
{
    return (head + 1) % size == tail;
}

/// @brief Check if there is any new data available in the buffer.
/// @return Whether or not new data is available.
bool RingBuffer::isDataAvailable()
{
    return !isEmpty;
}

/// @brief Get the last time that the buffer was marked as empty.
/// @return Time that buffer went empty.
chrono::time_point<chrono::high_resolution_clock> RingBuffer::getEmptyTime()
{
    return isEmptyTime;
}

/// @brief Get the current count of items in the buffer.
/// @return The current size of the buffer.
int RingBuffer::bufferSize()
{
    if (isEmpty)
    {
        return 0;
    }

    if (head == tail && !isEmpty)
    {
        return size;
    }

    if (head > tail)
    {
        return head - tail;
    }

    return size - (tail - head);
}

/// @brief Get the size that the buffer was initialized with.
/// @return The maximum size of the buffer.
int RingBuffer::maximumSize()
{
    return size;
}

/// @brief Print the current head and tail position of the buffer.
void RingBuffer::printStats()
{
    spdlog::info("Ring buffer head: {}, tail: {}", head, tail);
}

/// @brief Output the data that is currently in the buffer to the console.
void RingBuffer::printData()
{
    if (!isDataAvailable())
    {
        cout << "Buffer is empty!\n";

        return;
    }

    cout << "Content of buffer: ";

    int itemsInBuff = bufferSize();

    for (int i = 0; i < itemsInBuff; i++)
    {
        int element = buffer[(tail + i) % size];

        cout << element << (i < itemsInBuff - 1 ? ", " : "\n");
    }
}

/// @brief Check if the given received time is already present in the received times array:
/// @param receivedTime The time at which the current data was received.
/// @return Whether or not this time has been recorded already.
bool RingBuffer::receivedTimeSeen(const std::chrono::time_point<std::chrono::high_resolution_clock> &receivedTime)
{
    for (int i = 0; i < receivedTimes.size(); i++)
    {
        if (receivedTimes[i].second == receivedTime)
        {
            return true;
        }
    }

    return false;
}
