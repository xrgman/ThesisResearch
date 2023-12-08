#ifndef AUDIOCODEC_H
#define AUDIOCODEC_H

#include "main.h"

#define AUDIO_CODEC_SIZE 220500

struct AudioCodecResult {
    int senderId;
    double doa;
};

class AudioCodec
{
public:
    AudioCodec();

    void encode(int16_t *output, int outputSize, int senderId);

    AudioCodecResult decode(); //Input should be an array? Or we do it bit for bit, thats probably better

private:
    double volume;

    void generateChirp(int16_t *output, double startFrequency, double stopFrequency, double duration);
    void encodePreamble(int16_t *output, double startFrequency, double stopFrequency);

    int16_t applyKaiserWindow(int16_t value, int totalSize, int i, int beta);
};

#endif