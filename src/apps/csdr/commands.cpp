/*
Copyright (c) 2021-2023 Jakob Ketterl <jakob.ketterl@gmx.de>

This file is part of csdr.

csdr is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

csdr is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with csdr.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "commands.hpp"
#include "async.hpp"

#include "agc.hpp"
#include "fmdemod.hpp"
#include "amdemod.hpp"
#include "dcblock.hpp"
#include "converter.hpp"
#include "fft.hpp"
#include "logpower.hpp"
#include "logaveragepower.hpp"
#include "fftexchangesides.hpp"
#include "realpart.hpp"
#include "firdecimate.hpp"
#include "benchmark.hpp"
#include "fractionaldecimator.hpp"
#include "adpcm.hpp"
#include "limit.hpp"
#include "deemphasis.hpp"
#include "gain.hpp"
#include "dbpsk.hpp"
#include "varicode.hpp"
#include "timingrecovery.hpp"
#include "noise.hpp"
#include "phasedemod.hpp"
#include "rtty.hpp"
#include "baudot.hpp"
#include "mfrtty.hpp"
#include "sstv.hpp"
#include "fax.hpp"
#include "afc.hpp"
#include "cw.hpp"
#include "noisefilter.hpp"
#include "sitorb.hpp"
#include "ccir476.hpp"
#include "dsc.hpp"
#include "ccir493.hpp"
#include "navtex.hpp"

#include <iostream>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <cstdio>

using namespace Csdr;

template <typename T, typename U>
void Command::runModule(Module<T, U>* module) {
    auto buffer = new Ringbuffer<T>(bufferSize());
    module->setReader(new RingbufferReader<T>(buffer));
    auto writer = new StdoutWriter<U>();
    module->setWriter(writer);

    AsyncRunner* runner = nullptr;
    if (*get_parent()->get_option("--async")) {
        runner = new AsyncRunner(module);
    }

    fd_set read_fds;
    struct timeval tv = { .tv_sec = 10, .tv_usec = 0};
    int rc;
    size_t bytes_read;
    size_t read_over = 0;
    int nfds = fileno(stdin) + 1;

    FILE* fifo = nullptr;
    char* fifo_input = nullptr;
    if (!fifoName.empty()) {
        fifo = fopen(fifoName.c_str(), "r");
        if (fifo == nullptr) {
            std::cerr << "error opening fifo: " << strerror(errno) << "\n";
        } else {
            fcntl(fileno(fifo), F_SETFL, O_NONBLOCK);
            nfds = std::max(fileno(stdin), fileno(fifo)) + 1;
            fifo_input = (char*) malloc(1024);
        }
    }

    bool run = true;
    while (run) {
        FD_ZERO(&read_fds);
        FD_SET(fileno(stdin), &read_fds);
        if (fifo) FD_SET(fileno(fifo), &read_fds);
        tv.tv_sec = 10;
        tv.tv_usec = 0;

        rc = select(nfds, &read_fds, NULL, NULL, &tv);
        if (rc == -1) {
            std::cerr << "select() error: " << strerror(errno) << "\n";
            break;
        } else if (rc) {
            if (fifo && FD_ISSET(fileno(fifo), &read_fds)) {
                if (fgets(fifo_input, 1024, fifo) != NULL) {
                    processFifoData(std::string(fifo_input, strlen(fifo_input) - 1));
                } else {
                    std::cerr << "WARNING: fifo returned from select(), but no data.\n";
                }
            }
            if (FD_ISSET(fileno(stdin), &read_fds)) {
                // clamp so we don't overwrite the whole buffer in one go
                size_t writeable = std::min((size_t) 1024, buffer->writeable());
                // compensate for byte to element alignment
                writeable = (writeable * sizeof(T)) - read_over;
                bytes_read = read(fileno(stdin), ((char *) buffer->getWritePointer()) + read_over, writeable);
                if (bytes_read == 0) break;

                // advance but don't go into partially read elements
                buffer->advance((bytes_read + read_over) / sizeof(T));
                read_over = (bytes_read + read_over) % sizeof(T);

                // synchronous processing if we are not in async mode
                if (runner == nullptr) {
                    while (module->canProcess()) module->process();
                }
            }
        //} else {
            // no data, just timeout.
        }

        if (fifo && feof(fifo)) {
            std::cerr << "WARNING: fifo indicates EOF, terminating\n";
            run = false;
        }

    }

    if (fifo) {
        fclose(fifo);
        free(fifo_input);
    }

    delete buffer;
    delete runner;
}

template<typename T>
void Command::runSource(Source<T>* source) {
    auto writer = new StdoutWriter<T>();
    source->setWriter(writer);

    bool run = true;
    while (run) {
        sleep(1);
    }
}

CLI::Option* Command::addFifoOption() {
    return add_option("--fifo", fifoName, "Control fifo");
}

AgcCommand::AgcCommand(): Command("agc", "Automatic gain control") {
    add_set("-f,--format", format, {"s16", "float", "complex"}, "Data format", true);
    add_set("-p,--profile", profile, {"fast", "slow"}, "AGC profile", true);
    add_option("-a,--attack", attack, "AGC attack rate (slow: 0.1; fast: 0.01)");
    add_option("-d,--decay", decay, "AGC decay rate (slow: 0.0001; fast: 0.001)");
    add_option("-t,--hangtime", hangtime, "AGC hang time (slow: 600; fast: 200)");
    add_option("-m,--max", max_gain, "Maximum gain", true);
    add_option("-i,--initial", initial_gain, "Initial gain", true);
    add_option("-r,--reference", reference, "Reference level", true);
    callback( [this] () {
        if (format == "float") {
            runAgc<float>();
        } else if (format == "s16") {
            runAgc<short>();
        } else if (format == "complex") {
            runAgc<complex<float>>();
        } else {
            std::cerr << "invalid format: " << format << "\n";
        }
    });
}

template <typename T>
void AgcCommand::runAgc() {
    Agc<T>* agc = new Agc<T>();
    if (profile == "fast") {
        agc->setAttack(0.1);
        agc->setDecay(0.001);
        agc->setHangTime(200);
    } else if (profile == "slow") {
        agc->setAttack(0.01);
        agc->setDecay(0.0001);
        agc->setHangTime(600);
    }
    // override profile settings with explicit options on the command-line, if given
    if (attack != 0) agc->setAttack(attack);
    if (decay != 0) agc->setDecay(decay);
    if (hangtime != 0) agc->setHangTime(hangtime);
    agc->setMaxGain(max_gain);
    agc->setInitialGain(initial_gain);
    agc->setReference(reference);
    runModule(agc);
}

FmdemodCommand::FmdemodCommand(): Command("fmdemod", "FM demodulation") {
    callback( [this] () {
        runModule(new FmDemod());
    });
}

AmdemodCommand::AmdemodCommand(): Command("amdemod", "AM demodulation") {
    callback( [this] () {
        runModule(new AmDemod());
    });
}

DcBlockCommand::DcBlockCommand(): Command("dcblock", "DC block") {
    callback( [this] () {
        runModule(new DcBlock());
    });
}

ConvertCommand::ConvertCommand(): Command("convert", "Convert between stream formats") {
    add_set("-i,--informat", inFormat, {"s16", "float", "char"}, "Input data format", true);
    add_set("-o,--outformat", outFormat, {"s16", "float", "char"}, "Output data format", true);
    callback( [this] () {
        if (inFormat == outFormat) {
            std::cerr << "input and output format are identical, cannot convert\n";
            return;
        }
        if (inFormat == "s16") {
            if (outFormat == "float") {
                runModule(new Converter<short, float>());
            } else {
                std::cerr << "unable to handle output format \"" << outFormat << "\"\n";
            }
        } else if (inFormat == "float") {
            if (outFormat == "s16") {
                runModule(new Converter<float, short>());
            } else if (outFormat == "char") {
                runModule(new Converter<float, unsigned char>());
            } else {
                std::cerr << "unable to handle output format \"" << outFormat << "\"\n";
            }
        } else if (inFormat == "char") {
            if (outFormat == "float") {
                runModule(new Converter<unsigned char, float>());
            } else {
                std::cerr << "unable to handle output format \"" << outFormat << "\"\n";
            }
        } else {
            std::cerr << "unable to handle input format \"" << inFormat << "\"\n";
        }
    });
}

FftCommand::FftCommand(): Command("fft", "Fast Fourier transformation") {
    add_option("fft_size", fftSize, "FFT size")->required();
    add_option("every_n_samples", everyNSamples, "Run FFT every N samples")->required();
    add_set("-w,--window", window, {"boxcar", "blackman", "hamming"}, "Window function", true);
    callback( [this] () {
        if (!isPowerOf2(fftSize)) {
            std::cerr << "FFT size must be power of 2\n";
            return;
        }
        Window* w;
        if (window == "boxcar") {
            w = new BoxcarWindow();
        } else if (window == "blackman") {
            w = new BlackmanWindow();
        } else if (window == "hamming") {
            w = new HammingWindow();
        } else {
            std::cerr << "window type \"" << window << "\" not available\n";
            return;
        }

        runModule(new Fft(fftSize, everyNSamples, w));
    });
}

bool FftCommand::isPowerOf2(unsigned int size) {
    unsigned char bitcount = 0;
    for (int i = 0; i < 32; i++) {
        bitcount += (size >> i) & 1;
    }
    return bitcount == 1;
}

LogPowerCommand::LogPowerCommand(): Command("logpower", "Calculate dB power") {
    add_option("add_db", add_db, "Offset in dB", true);
    callback( [this] () {
        runModule(new LogPower(add_db));
    });
}

LogAveragePowerCommand::LogAveragePowerCommand(): Command("logaveragepower", "Calculate average dB power") {
    add_option("fft_size", fftSize, "Number of FFT bins")->required();
    add_option("avg_number", avgNumber, "Number of FFTs to average")->required();
    add_option("-a,--add", add_db, "Offset in dB", true);
    callback( [this] () {
        runModule(new LogAveragePower(fftSize, avgNumber, add_db));
    });
}

FftExchangeSidesCommand::FftExchangeSidesCommand(): Command("fftswap", "Switch FFT sides") {
    add_option("fft_size", fftSize, "Number of FFT bins")->required();
    callback( [this] () {
        runModule(new FftExchangeSides(fftSize));
    });
}

RealpartCommand::RealpartCommand(): Command("realpart", "Extract the real part of an IQ signal") {
    callback( [this] () {
        runModule(new Realpart());
    });
}

ShiftCommand::ShiftCommand(): Command("shift", "Shift a signal in the frequency domain") {
    add_option("rate", rate, "Amount of shift relative to the sampling rate");
    addFifoOption();
    callback( [this] () {
        //auto shift = new ShiftMath(rate);
        auto shift = new ShiftAddfast(rate);
        shiftModule = shift;
        runModule(shift);
    });
}

void ShiftCommand::processFifoData(std::string data) {
    shiftModule->setRate(std::stof(data));
}

FirDecimateCommand::FirDecimateCommand(): Command("firdecimate", "Decimate and filter") {
    add_option("decimation_factor", decimationFactor, "Decimation factor")->required();
    add_option("transition_bw", transitionBandwidth, "Transition bandwidth", true);
    add_option("-c,--cutoff", cutoffRate, "Cutoff rate", true);
    add_set("-w,--window", window, {"boxcar", "blackman", "hamming"}, "Window function", true);
    callback( [this] () {
        Window* w;
        if (window == "boxcar") {
            w = new BoxcarWindow();
        } else if (window == "blackman") {
            w = new BlackmanWindow();
        } else if (window == "hamming") {
            w = new HammingWindow();
        } else {
            std::cerr << "window type \"" << window << "\" not available\n";
            return;
        }
        runModule(new FirDecimate(decimationFactor, transitionBandwidth, w, cutoffRate));
    });
}

BenchmarkCommand::BenchmarkCommand(): Command("benchmark", "Perform internal benchmarks") {
    callback( [this] () {
        (new Benchmark())->run();
    });
}

FractionalDecimatorCommand::FractionalDecimatorCommand(): Command("fractionaldecimator", "Decimate in fractions") {
    add_set("-f,--format", format, {"float", "complex"}, "Format", true);
    add_option("decimation_rate", decimation_rate, "Decimation rate")->required();
    add_option("-n,--numpoly", num_poly_points, "Number of poly points", true);
    add_option("-t,--transition", transition, "Transition bandwidth for the prefilter", true);
    add_set("-w,--window", window, {"boxcar", "blackman", "hamming"}, "Window function for the prefilter", true);
    add_flag("-p,--prefilter", prefilter, "Apply filtering before decimation");
    callback( [this] () {
        if (format == "float") {
            runDecimator<float>();
        } else if (format == "complex") {
            runDecimator<complex<float>>();
        } else {
            std::cerr << "invalid format \"" << format << "\"\n";
        }
    });
}

template <typename T>
void FractionalDecimatorCommand::runDecimator() {
    LowPassFilter<T>* filter = nullptr;
    if (prefilter) {
        Window* w;
        if (window == "boxcar") {
            w = new BoxcarWindow();
        } else if (window == "blackman") {
            w = new BlackmanWindow();
        } else if (window == "hamming") {
            w = new HammingWindow();
        } else {
            std::cerr << "window type \"" << window << "\" not available\n";
            return;
        }
        filter = new LowPassFilter<T>(0.5 / (decimation_rate - transition), transition, w);
    }
    runModule(new FractionalDecimator<T>(decimation_rate, num_poly_points, filter));
}

AdpcmCommand::AdpcmCommand(): Command("adpcm", "ADPCM codec") {
    auto decodeFlag = add_flag("-d,--decode", decode, "Decode ADPCM data");
    auto encodeFlag = add_flag("-e,--encode", encode, "Encode into ADPCM data");
    // mutually exclusive
    encodeFlag->excludes(decodeFlag);
    decodeFlag->excludes(encodeFlag);
    add_flag("-s,--sync", sync, "Enable embedded sync frames");
    callback( [this] () {
        // default is encode
        if (!decode) {
            runModule(new AdpcmEncoder(sync));
        } else {
            runModule(new AdpcmDecoder());
        }
    });
}

FftAdpcmCommand::FftAdpcmCommand(): Command("fftadpcm", "Specialized ADPCM for FFT") {
    add_option("fft_size", fftSize, "Number of FFT bins")->required();
    callback( [this] () {
        runModule(new FftAdpcmEncoder(fftSize));
    });
}

LimitCommand::LimitCommand(): Command("limit", "Limit stream values to maximum amplitude") {
    add_option("max_amplitude", maxAmplitude, "Maximum amplitude", true);
    callback( [this] () {
        runModule(new Limit(maxAmplitude));
    });
}

PowerCommand::PowerCommand(): Command("power", "Measure power") {
    add_option("-o,--outfifo", outFifoName, "Control fifo")->required();
    add_option("length", length, "Number of samples to measure power over", true);
    add_option("decimation", decimation, "Decimate data when calculating power", true);
    add_option("report_every", reportInterval, "Reporting interval", true);
    callback( [this] () {
        int reportCounter = reportInterval;
        FILE* outFifo = fopen(outFifoName.c_str(), "w");
        if (outFifo == nullptr) {
            std::cerr << "error opening fifo: " << strerror(errno) << "\n";
            return;
        } else {
            fcntl(fileno(outFifo), F_SETFL, O_NONBLOCK);
        }
        runModule(new Power<complex<float>>(length, decimation, [this, &reportCounter, outFifo] (float power) {
            if (--reportCounter <= 0) {
                fprintf(outFifo, "%g\n", power);
                fflush(outFifo);
                reportCounter = reportInterval;
            }
        }));
        fclose(outFifo);
    });
}

SquelchCommand::SquelchCommand(): Command("squelch", "Measure power and apply squelch") {
    addFifoOption()->required();
    add_option("-o,--outfifo", outFifoName, "Control fifo")->required();
    add_option("length", length, "Number of samples to measure power over", true);
    add_option("decimation", decimation, "Decimate data when calculating power", true);
    add_option("hangLength", hangLength, "Number of samples to keep once power falls below threshold", true);
    add_option("flushLength", flushLength, "Number of samples to flush once squelch closes", true);
    add_option("report_every", reportInterval, "Reporting interval", true);
    callback( [this] () {
        int reportCounter = reportInterval;
        FILE* outFifo = fopen(outFifoName.c_str(), "w");
        if (outFifo == nullptr) {
            std::cerr << "error opening fifo: " << strerror(errno) << "\n";
            return;
        } else {
            fcntl(fileno(outFifo), F_SETFL, O_NONBLOCK);
        }
        squelch = new Squelch<complex<float>>(length, decimation, hangLength, flushLength, [this, &reportCounter, outFifo] (float power) {
            if (--reportCounter <= 0) {
                fprintf(outFifo, "%g\n", power);
                fflush(outFifo);
                reportCounter = reportInterval;
            }
        });
        runModule(squelch);
        fclose(outFifo);
    });
}

void SquelchCommand::processFifoData(std::string data) {
    squelch->setSquelch(std::stof(data));
}

SnrCommand::SnrCommand(): Command("snr", "Measure signal-to-noise ratio") {
    add_option("-o,--outfifo", outFifoName, "Control fifo")->required();
    add_option("length", length, "Number of samples to measure SNR over", true);
    add_option("fft_size", fftSize, "Size of the FFT being used", true);
    add_option("report_every", reportInterval, "Reporting interval", true);
    callback( [this] () {
        int reportCounter = reportInterval;
        FILE* outFifo = fopen(outFifoName.c_str(), "w");
        if (outFifo == nullptr) {
            std::cerr << "error opening fifo: " << strerror(errno) << "\n";
            return;
        } else {
            fcntl(fileno(outFifo), F_SETFL, O_NONBLOCK);
        }
        runModule(new Snr<complex<float>>(length, fftSize, [this, &reportCounter, outFifo] (float snr) {
            if (--reportCounter <= 0) {
                fprintf(outFifo, "%g\n", snr);
                fflush(outFifo);
                reportCounter = reportInterval;
            }
        }));
        fclose(outFifo);
    });
}

SnrSquelchCommand::SnrSquelchCommand(): Command("snrsquelch", "Measure signal-to-nouse ratio and apply squelch") {
    addFifoOption()->required();
    add_option("-o,--outfifo", outFifoName, "Control fifo")->required();
    add_option("length", length, "Number of samples to measure power over", true);
    add_option("fft_size", fftSize, "Size of the FFT being used", true);
    add_option("hangLength", hangLength, "Number of samples to keep once power falls below threshold", true);
    add_option("flushLength", flushLength, "Number of samples to flush once squelch closes", true);
    add_option("report_every", reportInterval, "Reporting interval", true);
    callback( [this] () {
        int reportCounter = reportInterval;
        FILE* outFifo = fopen(outFifoName.c_str(), "w");
        if (outFifo == nullptr) {
            std::cerr << "error opening fifo: " << strerror(errno) << "\n";
            return;
        } else {
            fcntl(fileno(outFifo), F_SETFL, O_NONBLOCK);
        }
        squelch = new SnrSquelch<complex<float>>(length, fftSize, hangLength, flushLength, [this, &reportCounter, outFifo] (float snr) {
            if (--reportCounter <= 0) {
                fprintf(outFifo, "%g\n", snr);
                fflush(outFifo);
                reportCounter = reportInterval;
            }
        });
        runModule(squelch);
        fclose(outFifo);
    });
}

void SnrSquelchCommand::processFifoData(std::string data) {
    squelch->setSquelch(std::stof(data));
}

DeemphasisCommand::DeemphasisCommand(): Command("deemphasis", "Deemphasis for FM applications") {
    auto wfmFlag = add_flag("-w,--wfm", "Wideband FM");
    auto nfmFlag = add_flag("-n,--nfm", "Narrowband FM");
    wfmFlag->excludes(nfmFlag);
    nfmFlag->excludes(wfmFlag);
    add_option("sample_rate", sampleRate, "Sample rate")->required();
    add_option("tau", tau, "Tau", true);
    callback( [this, wfmFlag] () {
        // default: nfm
        if (!(*wfmFlag)) {
            runModule(new NfmDeephasis(sampleRate));
        } else {
            runModule(new WfmDeemphasis(sampleRate, tau));
        }
    });
}

GainCommand::GainCommand(): Command("gain", "Apply fixed gain") {
    // there's no technical difference in float and complex gain, so we only implement float here
    add_option("gain", gain, "Gain factor")->required();
    callback( [this] () {
        runModule(new Gain<float>(gain));
    });
}

BandPassCommand::BandPassCommand(): Command("bandpass", "Bandpass filter") {
    addFifoOption();
    add_option("--low", lowcut, "Lower frequency");
    add_option("--high", highcut, "Higher Frequency");
    add_option("transition_bw", transition, "Transition bandwidth")->required();
    add_option("-w,--window", window, "Windowing function", true);
    add_flag("-f,--fft", use_fft, "Use FFT transformation filter");
    callback( [this] () {
        if (window == "boxcar") {
            windowObj = new BoxcarWindow();
        } else if (window == "blackman") {
            windowObj = new BlackmanWindow();
        } else if (window == "hamming") {
            windowObj = new HammingWindow();
        } else {
            std::cerr << "window type \"" << window << "\" not available\n";
            return;
        }
        if (use_fft) {
            auto filter = new FftBandPassFilter(lowcut, highcut, transition, windowObj);
            module = new FilterModule<complex<float>>(filter);
        } else {
            auto filter = new BandPassFilter<complex<float>>(lowcut, highcut, transition, windowObj);
            module = new FilterModule<complex<float>>(filter);
        }
        runModule(module);
    });
}

void BandPassCommand::processFifoData(std::string data) {
    std::stringstream ss(data);
    ss >> lowcut >> highcut;
    if (use_fft) {
        module->setFilter(new FftBandPassFilter(lowcut, highcut, transition, windowObj));
    } else {
        module->setFilter(new BandPassFilter<complex<float>>(lowcut, highcut, transition, windowObj));
    }
}

DBPskDecoderCommand::DBPskDecoderCommand(): Command("dbpskdecode", "Differential BPSK decoder") {
    callback( [this] () {
        runModule(new DBPskDecoder());
    });
}

VaricodeDecoderCommand::VaricodeDecoderCommand(): Command("varicodedecode", "Decode PSK31 varicode") {
    callback( [this] () {
        runModule(new VaricodeDecoder());
    });
}

TimingRecoveryCommand::TimingRecoveryCommand(): Command("timingrecovery", "Timing recovery") {
    add_set("-f,--format", format, {"float", "complex"}, "Data format", true);
    add_set("-a,--algorithm", algorithm, {"gardner", "earlylate"}, "Algorithm to be used", true);
    add_option("decimation", decimation, "Decimation (samples per symbol)")->required();
    add_option("loop_gain", loop_gain, "Loop gain", true);
    add_option("max_error", max_error, "Max allowed error", true);
    callback( [this] () {
        if (algorithm == "gardner") {
            if (format == "float") {
                runModule(new GardnerTimingRecovery<float>(decimation, loop_gain, max_error));
            } else if (format == "complex") {
                runModule(new GardnerTimingRecovery<complex<float>>(decimation, loop_gain, max_error));
            } else {
                std::cerr << "Invalid format: \"" << format << "\"\n";
            }
        } else if (algorithm == "earlylate") {
            if (format == "float") {
                runModule(new EarlyLateTimingRecovery<float>(decimation, loop_gain, max_error));
            } else if (format == "complex") {
                runModule(new EarlyLateTimingRecovery<complex<float>>(decimation, loop_gain, max_error));
            } else {
                std::cerr << "Invalid format: \"" << format << "\"\n";
            }
        } else {
            std::cerr << "Invalid algorithm: \"" << algorithm << "\"\n";
        }
    });
}

NoiseCommand::NoiseCommand(): Command("noise", "Noise generator") {
    callback([this] () {
        runSource(new NoiseSource<complex<float>>());
    });
}

Phasedemodcommand::Phasedemodcommand(): Command("phasedemod", "Phase demodulation") {
    callback([this] () {
        runModule(new PhaseDemod());
    });
}

RttyDecodeCommand::RttyDecodeCommand(): Command("rttydecode", "RTTY decoder") {
    add_flag("-i,--invert", invert, "Inverse operation (swap MARK/SPACE)");
    callback([this] () {
        runModule(new RttyDecoder(invert));
    });
}

BaudotDecodeCommand::BaudotDecodeCommand(): Command("baudotdecode", "Baudot decoder") {
    callback([this] () {
        runModule(new BaudotDecoder());
    });
}

LowpassCommand::LowpassCommand(): Command("lowpass", "Lowpass FIR filter") {
    add_set("-f,--format", format, {"float", "complex"}, "Data format", true);
    add_option("cutoff", cutoffRate, "Cutoff rate")->required();
    add_option("transition_bw", transitionBandwidth, "Transition bandwidth", true);
    add_set("-w,--window", window, {"boxcar", "blackman", "hamming"}, "Window function", true);
    callback([this] () {
        Window* w;
        if (window == "boxcar") {
            w = new BoxcarWindow();
        } else if (window == "blackman") {
            w = new BlackmanWindow();
        } else if (window == "hamming") {
            w = new HammingWindow();
        } else {
            std::cerr << "window type \"" << window << "\" not available\n";
            return;
        }
        if (format == "float") {
            runModule(new FilterModule<float>(new LowPassFilter<float>(cutoffRate, transitionBandwidth, w)));
        } else if (format == "complex") {
            runModule(new FilterModule<complex<float>>(new LowPassFilter<complex<float>>(cutoffRate, transitionBandwidth, w)));
        } else {
            std::cerr << "invalid format: " << format << "\n";
        }
    });
}

CwDecoderCommand::CwDecoderCommand(): Command("cwdecode", "CW decoder") {
    add_option("sample_rate", sampleRate, "Sample rate")->required();
    add_flag("-c,--cw", showCw, "Show CW code (dits and dahs)");
    callback( [this] () {
        runModule(new CwDecoder<float>(sampleRate, showCw));
    });
}

MFRttyDecoderCommand::MFRttyDecoderCommand(): Command("mfrttydecode", "RTTY decoder") {
    add_option("sample_rate", sampleRate, "Sample rate")->required();
    add_option("freq", targetFreq, "Frequency base");
    add_option("shift", targetWidth, "Frequency shift");
    add_option("baud_rate", baudRate, "Baud rate");
    add_option("reverse", reverse, "Reverse space and mark");
    callback( [this] () {
        runModule(new MFRttyDecoder<float>(sampleRate, targetFreq, targetWidth, baudRate, reverse));
    });
}

SstvDecoderCommand::SstvDecoderCommand(): Command("sstvdecode", "SSTV decoder") {
    add_option("sample_rate", sampleRate, "Sample rate")->required();
    callback( [this] () {
        runModule(new SstvDecoder<float>(sampleRate));
    });
}

FaxDecoderCommand::FaxDecoderCommand(): Command("faxdecode", "FAX decoder") {
    add_option("sample_rate", sampleRate, "Sample rate")->required();
    add_option("lpm", lpm, "Lines per minute");
    add_option("max_lines", maxLines, "Maximum number of lines");
    add_option("am", am, "Use AM modulation");
    add_option("color", color, "Color fax");
    add_option("sync", sync, "Sync scanlines");

    unsigned int options =
        (am?    FaxDecoder<float>::OPT_AM    : 0) |
        (sync?  FaxDecoder<float>::OPT_SYNC  : 0) |
        (color? FaxDecoder<float>::OPT_COLOR : 0);

    callback( [this, options] () {
        runModule(new FaxDecoder<float>(sampleRate, lpm, maxLines, options));
    });
}

ReduceNoiseCommand::ReduceNoiseCommand(): Command("reducenoise", "Reduce noise") {
    addFifoOption();
    add_option("-f,--fft_size", fftSize, "Number of FFT bins");
    add_option("-w,--wnd_size", wndSize, "Filter window size");
    add_option("-l,--latency", latency, "Latency");
    add_option("-t,--threshold", dBthreshold, "Suppression threshold in dB");
    callback( [this] () {
        auto filter = new AFNoiseFilter(fftSize, wndSize, latency);
        module = new FilterModule<float>(filter);
        filter->setThreshold(dBthreshold);
        runModule(module);
    });
}

AfcCommand::AfcCommand(): Command("afc", "Automatic frequency control") {
    add_option("update_period", updatePeriod, "Update period (>= sample_period)");
    add_option("sample_period", samplePeriod, "Sample period (>= 1)");

    callback( [this] () {
        runModule(new Afc(updatePeriod, samplePeriod));
    });
}

SitorBDecodeCommand::SitorBDecodeCommand(): Command("sitorbdecode", "SITOR-B decoder") {
    add_option("-e,--errors_allowed", errorsAllowed, "Number of errors allowed");
    add_flag("-i,--invert", invert, "Inverse operation (swap MARK/SPACE)");
    callback([this] () {
        runModule(new SitorBDecoder(errorsAllowed, invert));
    });
}

Ccir476DecodeCommand::Ccir476DecodeCommand(): Command("ccir476decode", "CCIR476 decoder") {
    callback([this] () {
        runModule(new Ccir476Decoder());
    });
}

DscDecodeCommand::DscDecodeCommand(): Command("dscdecode", "DSC decoder") {
    callback([this] () {
        runModule(new DscDecoder());
    });
}

Ccir493DecodeCommand::Ccir493DecodeCommand(): Command("ccir493decode", "CCIR493 decoder") {
    add_option("-e,--errors_allowed", errorsAllowed, "Number of errors allowed");
    add_flag("-i,--invert", invert, "Inverse operation (swap MARK/SPACE)");
    callback([this] () {
        runModule(new Ccir493Decoder(errorsAllowed, invert));
    });
}

NavtexDecodeCommand::NavtexDecodeCommand(): Command("navtexdecode", "NAVTEX decoder") {
    callback([this] () {
        runModule(new NavtexDecoder());
    });
}
