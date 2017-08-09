JPSTH ANALYSIS, 2008

written by: John Haitas, Jeremiah Cohen, and Jeff Schall
________________________________________________________________________

INPUT DATA FORMAT:

Spike data should be in the form of a matrix of spike times (spikeData), wherein each row denotes a separate trial.
_________________________________________________________________________

HOW TO RUN:

-Use the command:
 
	alignedSpikeData = alignTimeStamps(spikeData, alignEventVector) 

for each signal you wish to compare. This function aligns spike data to the stimulus of interest. The variable alignEventVector should be in the form of a column vector listing the event time for each trial.

-Use the command:

	timeStamps = trimTimeStamps(alignedSpikeData, timeWindow) 

for each signal you wish to compare. This function removes spike times outside the window of interest (timeWindow). The variable timeWindow should be in the form of [start finish], where both times are relative to the align event.

-Use the command:

	spikes = spikeCounts(timeStamps, timeWindow, binWidth) 

for each signal you wish to compare. This function counts the number of spikes in each bin for each trial.

- Run jpsth(spikeCounts_signal1, spikeCounts_signal2, coincidenceHistogramWidth). The output is a struct, which is described in the next section. 

________________________________________________________________________
OUTPUT STRUCTURE
________________________________________________________________________

jpsthData 		
		-> superstructure for jpsth output data

jpsthData.psth_1 
		-> peristimulus time histogram for signal 1

jpsthData.psth_2 
		-> peristimulus time histogram for signal 2

jpsthData.normalizedJPSTH
		-> jpsth matrix normalized to the standard deviation

jpsthData.xcorrHist
		-> crosscorrelogram based on Aertsen et al's equations

jpsthData.pstch
		-> coincidence histogram based on Aertsen et al's equations

jpsthData.covariogram
		-> covariogram based on Brody's equation (qualitatively same as crosscorrelogram)

jpsthData.sigLow
		-> negative covariogram significance limit based on Brody's equation

jpsthData.sigHigh
		-> positive covariogram significance limit based on Brody's equation

jpsthData.sigPeakEndpoints 
		-> start and end indices for continuous span of covariogram values exceeding sigHigh

jpsthData.sigTroughEndpoints 
		-> start and end indices for continuous span of covariogram values less than sigLow
