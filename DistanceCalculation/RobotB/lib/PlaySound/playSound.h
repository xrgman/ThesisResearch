#ifndef PLAYSOUND_H
#define PLAYSOUND_H

#include "mbed.h"

Timer outputTimer;
bool playingSound = false;

/// @brief Play a beep over the speaker.
/// @param speaker PWM object of the speaker output.
/// @return Whether or not playing is done.
bool playBeep(PwmOut &speaker)
{
  if (!playingSound)
  {
    //speaker.write(0.5f);

    outputTimer.reset();
    outputTimer.start();

    playingSound = true;

    return false;
  }

  auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(outputTimer.elapsed_time()).count();

  if (elapsedTime >= 500)
  {
    //speaker.write(0.0f);

    outputTimer.stop();
    playingSound = false;

    return true;
  }

  return false;
}


#endif