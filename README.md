# ThesisResearch

## TO-DO implementation:

<ol>
  <li><s>Play sound over speaker (WAV file)</s></li>
  <li><s>Record audio from microphones (WAV file)</s></li>
  <li><s>Show plot of frequency spectrum</s></li>
  <li><s>Correct the microphone order</s></li>
  <li><s>Create the possibility to load in map</s></li>
  <li><s>Create basic particle filter, by plotting dots uniformly</s></li>
  <li><s>Update particle filter on movement</s></li>
  <li><s>Show current cell based on particles</s></li>
  <li><s>Encode message data in the form of audio</s></li>
  <li><s>Filter out messages from own speaker</s></li>
  <li><s>Create message protocol (chirp) between robots via audio</s></li>
  <li><s>Determine DOA of received chirp</s></li>
  <li><s>Add the I've seen a wall message type functionality</s></li>
  <ol type="1">
      <li><s>Sends distance to wall.</s></li>
      <li><s>Sends orientation of wall with respect to north.</s></li>
  </ol>
  <li><s>Add the I'm in this cell message type functionality</s></li>
   <ol type="1">
      <li><s>Sends cell the robot is in.</s></li>
  </ol>
  <li>Check that we do not process message from own ID</li>
  <li>Add the ability to process wall message and update particle filter based on it</li>
  <li>Add the ability to process cell message and update particle filter based on it</li>
  <li><s>Cache distances between cells to a .json file to speed up particle filter</s></li>
  <li>Save the localization table to file after generating it</li>
  <li>Add the ability to load a localization table from another robot and fuse it with own tables. If none are there, set it as the table for that robot.</li>
  <li>Increase performance PF with high number of cells</li>
  <li>Link movement from robot to movement in particle filter</li> 
</ol>

## TO-DO evaluation
<ol>
  <li>Evaluate DOA performance and based on findings maybe update it even more</li>
  <li>Compare localization accaracy compared to normal PF operation</li>
  <li>Evaluate distance performance in terms of BER and against number of samples used in preamble, per bit</li>
  <li>Compare communication protocol performance with other works</li>
   <ol type="1">
      <li>Compare SNR performance for both white gaussian noise and babble noise.</li>
      <li>Compare bitrate.</li>
      <li>Compare distance performance (without extra noise, but with recordings), in terms of BER.</li>
      <li>Compare multi robot per distance performance, in terms of BER.</li>
      <li>Compare localization accuracy, in terms of wrong cell?.</li>
    </ol>
</ol>

## Algorithm:
<ol>
  <li>Robot drives to wall on the right side of its room (with respect to north of map)</li>
    <ol type="1">
      <li>Particles are being updated during driving, just like normal PF operation.</li>
      <li>Walls can be detected by continues obstacle in the place we would expect a wall, objects in the room won't have this property. Advanced wall detection is out of scope?</li>
      <li>When wall is reached & detected, eleminate all particles that are not within x cm from a right wall.</li>
    </ol>
  <li>Send message to other robots and receive messages from other robots.</li>
  <ol type="1">
      <li>Based on the DOA of the message more particles should be removed and we should check for convergence.</li>
      <li>If no convergence yet, drive to top or bottom wall and repeat the progress.</li>
      <li>If convergence, send out I'm in this cell message to other robots.</li>
  </ol>
  <li>Receive I'm in this cell message from other robot.</li>
  <ol type="1">
      <li>Based on the DOA of the message, a lot of cells can be elimanted and convergence should be achieved.</li>
  </ol>
  <li>Upon driving keep updating particles and processing messages and sending out I'm in this cell messages.</li>
</ol>

## Assumptions:
<ol>
  <li>IMU implementation is out of scope and we assume we can just get the current orientation of the robot.</li>
  <li>Driving functionality is implemented, but during the evaluation robot will be picked up and moved for more accurate evaluation and to make sure driving won't be a bottleneck.</li>
  <li>Wall detection is assumed to be possible (mention a ref to paper in thesis), advanced will not be implemented as we will simply tell the robot the wall is x cm away for evaluation purposes.</li>
  <li>Angle of wall and distance to wall have 3 decimals.</li>
</ol>


## Installation libraries:
git clone https://github.com/PortAudio/portaudio.git
./configure && make
sudo make install
sudo ln -s /usr/local/lib/libportaudio.so /usr/lib/arm-linux-gnueabihf/libportaudio.so

sudo ln -s /usr/local/lib/libportaudio.so /usr/lib/aarch64-linux-gnu/libportaudio.so

sudo apt-get install libspdlog-dev libboost-all-dev libsdl2-ttf-dev libsdl2-dev

Use alsamixer on rpi to increase volume to <b>81</b> of speaker.

### Stop tracking changes to config:
git update-index --assume-unchanged src/config.json

### Start tracking changes to config:
git update-index --no-assume-unchanged src/config.json