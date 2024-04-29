#ifndef SOUNDDETECTION_H
#define SOUNDDETECTION_H

#include "mbed.h"
#include "../Util/util.h"

// Function definitions:
void resetBeepDetection();

uint16_t previousValue = 0;

/// @brief Detect if there is a beep in the room, based on the current and previous sample.
/// @param sample New sample. 
/// @return Whether or not a beep was detected.
bool detectBeep(uint16_t sample)
{
    uint16_t peakToPeak = 0;

    // Calculating peak to peak value:
    if (previousValue > 0)
    {
        peakToPeak = abs(previousValue - sample);
    }

    // Saving current sample:
    previousValue = sample;

    if (peakToPeak > 400) //450
    {
        return true;
    }

    return false;
}

/// @brief Resetting all values used in beep detection:
void resetBeepDetection()
{
    previousValue = 0;
}

#endif