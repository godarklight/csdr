/*
This file is part of libcsdr.

	Copyright (c) Andras Retzler, HA7ILM <randras@sdr.hu>
	Copyright (c) Warren Pratt, NR0V <warren@wpratt.com>
    Copyright (c) Jakob Ketterl, DD5JFK <jakob.ketterl@gmx.de>
	Copyright 2006,2010,2012 Free Software Foundation, Inc.

    libcsdr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libcsdr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libcsdr.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "agc.hpp"
#include "complex.hpp"

#include <cmath>
#include <climits>
#include <algorithm>
#include <iostream>

using namespace Csdr;

template <typename T>
void Agc<T>::process(T* input, T* output, size_t work_size) {
    float input_abs, envelope, target_gain;

    //Sets the amount of weight target_gain affects each sample.
    //Make sure that you can get close the the target by lookahead_size samples.
    float beta = 0.02;

    for (int i = 0; i < work_size; i++) {
        //Grab the envelope of the signal to work out our target gain
        input_abs = this->abs(input[i]);
        envelope = this->getEnvelope(input_abs);
        double new_target_gain = reference / envelope;

        //If we have detected an event that would cause us to go over our reference, let's pause the agc here and force set target_gain.
        if ((input_abs * target_gain) > reference)
        {
            hang_samples = hang_time + lookahead_size;

            //We only want to allow decreases in gain here, this would otherwise be incorrect if we have a smaller impulse after an initial impulse.
            if (new_target_gain < target_gain)
            {
                target_gain = new_target_gain;
            }

        }

        //If we have detected an event that causes us to go over reference, pause the AGC level as it is directly set above.
        //If we would go over reference, the counter is reset along with the target gain above.
        if (hang_samples > 0)
        {
            hang_samples--;
        }
        else
        {
            //Set the target gain so that our signal is brought up to the reference level.
            target_gain = target_gain * (1.0 - beta) + new_target_gain * beta;
        }


        //Clamp gain to max_gain and 0
        if (target_gain > max_gain) target_gain = max_gain;
        if (target_gain < 0) target_gain = 0;

        //Lowpass filter the actual gain.
        gain = gain * (1.0 - beta) + (target_gain * beta);

        //std::cerr << "gain: " << gain << ", target_gain: " << target_gain << ", envelope: " << envelope << std::endl;

        //Because we have added delay, our first samples are in last_samples, and then the current input array is used.
        if (i < lookahead_size)
        {
            output[i] = scale(last_samples[i]);
        }
        else
        {
            output[i] = scale(input[i - lookahead_size]);
        }
    }
    //We need to save the last lookahead_size samples for the next call.
    for (int i = 0; i < lookahead_size; i++)
    {
        last_samples[i] = input[(work_size - lookahead_size) + i];
    }
}

template <typename T>
float Agc<T>::getEnvelope(float in)
{
    if (in > env_detect)
    {
        env_detect = in;
    }
    else
    {
        env_detect = env_detect * (1.0 - decay_rate);
    }
    return env_detect;
}


template <>
float Agc<short>::abs(short in) {
    return std::fabs((float) in) / SHRT_MAX;
}

template <>
bool Agc<short>::isZero(short in) {
    return in == 0;
}

template <>
short Agc<short>::scale(short in) {
    float val = gain * in;
    if (val >= SHRT_MAX) return SHRT_MAX;
    if (val <= SHRT_MIN) return SHRT_MIN;
    return (short) val;
}

template <>
float Agc<float>::abs(float in) {
    return std::fabs(in);
}

template <>
bool Agc<float>::isZero(float in) {
    return in == 0.0f;
}

template <>
float Agc<float>::scale(float in) {
    float val = in * gain;
    if (val > 1.0f) return 1.0f;
    if (val < -1.0f) return -1.0f;
    return val;
}

template <>
float Agc<complex<float>>::abs(complex<float> in) {
    return std::abs(in);
}

template <>
bool Agc<complex<float>>::isZero(complex<float> in) {
    return in == complex<float>(0, 0);
}

template <>
complex<float> Agc<complex<float>>::scale(complex<float> in) {
    complex<float> val = in * gain;
    //Don't clip IQ individually, clip the magnitude of the complex number. std::abs returns magnitude for std::complex.
    //This form of clipping sounds much better than clipping a real valued signal.
    float mag = abs(val);
    if (mag > 1.0)
    {
        val = val / mag;
    }
    return val;
}

template <typename T>
void Agc<T>::setReference(float reference) {
    this->reference = reference;
}

template <typename T>
void Agc<T>::setAttack(float attack_rate) {
    this->attack_rate = attack_rate;
}

template <typename T>
void Agc<T>::setDecay(float decay_rate) {
    this->decay_rate = decay_rate;
}

template <typename T>
void Agc<T>::setMaxGain(float max_gain) {
    this->max_gain = max_gain;
}

template <typename T>
void Agc<T>::setInitialGain(float initial_gain) {
    gain = initial_gain;
}

template <typename T>
void Agc<T>::setHangTime(unsigned long int hang_time) {
    this->hang_time = hang_time;
}

namespace Csdr {
    template class Agc<short>;
    template class Agc<float>;
    template class Agc<complex<float>>;
}