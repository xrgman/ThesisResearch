/* mbed Microcontroller Library
 * Copyright (c) 2019 ARM Limited
 * SPDX-License-Identifier: Apache-2.0
 */
#include "mbed.h"
#include "playSound.h"
#include "soundDetection.h"

using namespace std::chrono;

#define SPEAKER_PIN PA_5
#define SPEAKER_SELECT_PIN PC_4
#define MIC_PIN PC_3

#define SPEAKER_FREQUENCY 2000 // In Hz
#define SAMPLING_FREQ 44100 // In Hz

#define NR_OF_MESSAGES 99

// Function definitions:
void handleChar(char c);
void reset();

// Object definitions (3mbs is limit of uart chip on board):
UnbufferedSerial pc(USBTX, USBRX, 3000000);

AnalogIn microphone(MIC_PIN);
DigitalOut speakerSelectPin(SPEAKER_SELECT_PIN);
PwmOut speaker(SPEAKER_PIN);

// Field definitions:
bool newDataAvailable = false;

bool messageSend = false;
bool playSound = false;

// For message sending:
int nrOfMessagesSend = 0;
const int messageTime = SAMPLING_FREQ * 2; //2 seconds


int newDataCnt = 0;
int sampleCounter = 1;

bool found = false;

void microphoneEventHandler()
{
    if (messageSend)
    {
        newDataAvailable = true;
        newDataCnt++;

        sampleCounter++;
    }
}

int main()
{
    printf("********************\r\n");
    printf("Robot A started\r\n");
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

    // Serial events:
    pc.attach([]()
              {
        char command;
        pc.read(&command, 1); 
        handleChar(command); },
              UnbufferedSerial::RxIrq);

    while (true)
    {
        if (playSound)
        {
            if (playBeep(speaker))
            {
                found = false;
                playSound = false;


                sampleCounter = 1;
                newDataCnt = 0;
                messageSend = true;
            }
        }

        // Sending message:
        if (messageSend)
        {
            //Checking if new message should be send:
            if (sampleCounter % messageTime == 0)
            {
                messageSend = false;
                newDataAvailable = false;

                if (nrOfMessagesSend < NR_OF_MESSAGES)
                {
                    playSound = true;

                    nrOfMessagesSend++;
                }
                else
                {
                    nrOfMessagesSend = 0;

                    printf("Done!\r\n");
                }
            }

            // Handling microphone events:
            if (newDataAvailable)
            {
                uint16_t sample = microphone.read_u16();

                if (detectBeep(sample) && !found && sampleCounter > 20000)
                {
                    printf("%d\r\n", newDataCnt);

                    newDataCnt = 0;
                    found = true;
                }

                newDataAvailable = false;
            }
        }
    }
}

/// @brief Simple serial charachter processor.
/// @param c Characther received by serial event.
void handleChar(char c)
{
    switch (c)
    {
    case 's':
        messageSend = false;
        playSound = true;
        break;
    case 'r':
        reset();
        break;
    case 'p':
        reset();
        break;
    default:
        // pc.write("Input char not supported!\n", 26);
        break;
    }
}

/// @brief Reset all values back to their defaults.
void reset()
{
    messageSend = false;

    nrOfMessagesSend = 0;

    resetBeepDetection();
}
