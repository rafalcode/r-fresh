# From https://reintech.io/blog/biomedical-signal-processing-with-r-tutorial
# it's a bit poor

library(Cairo)
library(signal)
library(biosignalEMG)

# Let’s take an electrocardiogram (ECG) signal as an example. Here we will load the ECG signal, and then visualize it:

data(ecg)

CairoPNG("ecg00.png", 800, 800)
plot(ecg)
dev.off()

# This will display the ECG signal. Analyzing the signal might involve noise removal, filtering, and feature extraction. Here is an example of how to apply a low-pass filter to the ECG signal:
filtered_ecg <- butterworth(ecg, type="low", cutoff=100, order=4)
CairoPNG("ecg01.png", 800, 800)
plot(filtered_ecg)
dev.off()
