#include "main.h"

#include "gnuplot-iostream.h"

using namespace std;

void graphSineWave5()
{
    // Parameters for the sine wave
    const double frequency = 5000.0; // 5 kHz
    const double duration = 0.01;    // 0.01 seconds
    const int numPoints = 1000;      // Number of points to generate

    // Generate sine wave data
    vector<pair<double, double>> data;

    for (int i = 0; i < numPoints; ++i)
    {
        double time = i * duration / numPoints;
        double amplitude = sin(2.0 * M_PI * frequency * time);
        data.push_back(std::make_pair(time, amplitude));
    }

    // Set up Gnuplot
    Gnuplot gp;

    // Plot the sine wave
    gp << "set title 'Sine Wave (5 kHz)'\n";
    gp << "set xlabel 'Time (s)'\n";
    gp << "set ylabel 'Amplitude'\n";
    gp << "plot '-' with lines title 'Sine Wave'\n";
    gp.send1d(data);
}