#pragma once

class Stim {
public:
    double start = 0; // in ms
    double end = 0; // in ms
    double freq = 0; // in Hz
    double width = 0; // in ms

	/**
     * @brief Whether the stimulation is on at the given time.
     */
    bool pulse(double time);
};