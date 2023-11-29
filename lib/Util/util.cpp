#include "util.h"
#include <numeric>

using namespace std;

double calculateAverage(uint16_t *data, uint16_t size) {
    double sum = accumulate(data, data + size, 0.0);

    return sum / size;
}