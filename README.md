# SASP - Sound Analysis, Synthesis and Processing

This course covers a wide range of advanced topics of digital signal processing for the analysis of audio signals. Using signal processing concepts and tools such as windowing, overlap-and-add, Short-Time Fourier Transforms, uniform and iterated filterbanks, polyphase processing, time-space processing and statistical signal processing, we approach a variety of problems that are typical of sound analysis applications such as feature extraction and sound annotation, digital audio restoration, acoustic beamforming, source separation and localization. The course also addresses topics on digital signal processing for the synthesis/generation of timbres and sounds; for sound spatialization/reverberation; and for a high-end rendering of audio signals. Using signal processing concepts and tools such as Short-Time Fourier Transforms, uniform and iterated filterbanks, polyphase processing and space-time processing, we discuss a variety of applications such as musical sound synthesis, digital audio enhancement and improvement, multi-channel processing, acoustic beamforming and projection using speaker arrays, active acoustic conditioning, wavefield synthesis and holo-acoustics

HW1: The task is to perform short-time LPC analysis on both speech and music signals and feed the residual error signal of each music
segment to the all-pole filter computed from every speech segment. This way, the “excitation” pertaining to the music
signal (which contains pitch information) will be “shaped” by the time-varying resonant peaks characterizing the speech production
mechanism (e.g., the formants.).

HW2: Implementation of an acoustic source separation by applying binary masking to the short-time Fourier transform (STFT) of a mixture signal.

HW3: Implementation of a computationally efficient Leslie speaker emulation.

HW4: Modelling a three-way crossover network in the Wave Digital (WD) domain starting from a reference analog circuit.
