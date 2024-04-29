/* mbed Microcontroller Library
 * Copyright (c) 2019 ARM Limited
 * SPDX-License-Identifier: Apache-2.0
 */

#include "mbed.h"
#include "playSound.h"
#include "soundDetection.h"
#include "util.h"

using namespace std::chrono;

#define SPEAKER_PIN PA_5
#define SPEAKER_SELECT_PIN PC_4
#define MIC_PIN PC_3

#define SPEAKER_FREQUENCY 2000 // In Hz
#define SAMPLING_FREQ 44100    // In Hz

void reset();

// Object definitions:
UnbufferedSerial pc(USBTX, USBRX, 3000000);

AnalogIn microphone(MIC_PIN);
DigitalOut speakerSelectPin(SPEAKER_SELECT_PIN);
PwmOut speaker(PA_5);

// Field definitions:
bool newDataAvailable = false;

bool record = true;
bool playSound = false;

Timer processTimer;
unsigned long startTime, stopTime, peakDetectedTime;

int numberOfSamplesSeen = 0;

bool printSampleTime = false;

void microphoneEventHandler()
{
    if (record)
    {
        newDataAvailable = true;
    }
}

int main()
{
    printf("********************\r\n");
    printf("Robot B started\r\n");
    printf("********************\r\n");
    printf("\r\n");

    // Enabling amplifier on the board:
    speakerSelectPin = 1;

    // Set the PWM frequency to 10 kHz
    speaker.period(1.0f / SPEAKER_FREQUENCY); // Period is the inverse of the frequency

    // Set reference voltage for the microphone:
    microphone.set_reference_voltage(3.3f);

    Ticker ticker;
    ticker.attach(&microphoneEventHandler, microseconds(1000000) / SAMPLING_FREQ);

    reset();

    while (true)
    {
        if (playSound)
        {
            if (playBeep(speaker))
            {
                found = false;
                playSound = false;

                // Waiting half a second before starting listening again:
                ThisThread::sleep_for(chrono::milliseconds(500));

                reset();

                record = true;
            }
            else
            {
                continue;
            }
        }

        // Handling microphone events:
        if (newDataAvailable)
        {
            uint16_t sample = microphone.read_u16();
            ;
            if (numberOfSamplesSeen > 20)
            {
                if (detectBeep(sample))
                {
                    record = false;

                    // Waiting a second before sending:
                    ThisThread::sleep_for(chrono::milliseconds(1000));

                    playSound = true;
                }
            }

            newDataAvailable = false;
            numberOfSamplesSeen++;
        }
    }
}

/// @brief Reset all values back to their defaults.
void reset()
{
    numberOfSamplesSeen = 0;

    resetBeepDetection();
}
