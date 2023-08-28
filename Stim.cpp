#include "Stim.h"

bool Stim::pulse(double time) {
    if (time >= end || time < start)
        return false;
    double s = time - start; // time elapsed from the start
    int n = (int) (s * freq * 1e-3); // number of pulses that have happened or started
    if (freq) 
        s = time - start - n / freq * 1e3; // time elapsed in the current period.
    return (s < width);
}