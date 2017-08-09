function unitTest()
	
	windowMin = -10;
	windowMax = 30;
	window = [windowMin windowMax];
	
	binWidth = 10;
	
	signal1 = [-10 0 20; -25 -10 20; -45 30 45];
	signal2 = [-11 0 20; -22 -10 23; -45 30 43];
	
	expectedSpikes1 = [1 0 1 0; 0 0 1 0; 0 0 0 0];
	expectedSpikes2 = [1 0 1 0; 0 0 0 1; 0 0 0 0];
	
	expectedPSTH1 = [1/3 0 2/3 0];
	expectedPSTH2 = [1/3 0 1/3 1/3];	 
	
	% Run unit tests
	try
		% Trim time stamps
		[trimmedSignal1] = test_trimTimeStamps(signal1,window);
		[trimmedSignal2] = test_trimTimeStamps(signal2,window);

		% Test spikeCounts function
		[spikesSignal1] = test_spikeCounts(trimmedSignal1,window,binWidth,expectedSpikes1);
		[spikesSignal2] = test_spikeCounts(trimmedSignal2,window,binWidth,expectedSpikes2);

		% Test psth function
		[psth1] = test_psth(expectedSpikes1,expectedPSTH1);
		[psth2] = test_psth(expectedSpikes2,expectedPSTH2);

		% Test significantSpan function
		vector = [.1 .1 .1 .2 .2 .2 .1 .1 .1];
		sig = .15;
		thisSpan = [4 6];
		test_significantSpan(vector, sig, thisSpan)

		vector = [.1 .1 .2 .1 .1 .1 .1 .1 .1];
		sig = .15;
		thisSpan = [3 3];
		test_significantSpan(vector, sig, thisSpan)


		vector = [.1 .1 .1 .1 .1 .1 .1 .1 .1];
		sig = .15;
		thisSpan = [];
		test_significantSpan(vector, sig, thisSpan)
		  
		fprintf('%s\n','Unit Test SUCCESS')
	catch Exception
		beep
		fprintf('%s\n','Unit Test FAILURE');
		fprintf('%s\n',Exception.message);
		unitTestFailure = MException('VerifyOutput:OutOfBounds',...
						'Results are outside the allowable limits');
		throw(unitTestFailure);
	end
end

function [trimmedSignal] = test_trimTimeStamps(signal,window)
	trimmedSignal = trimTimeStamps(signal,window);
	try
  		assert(all(all(((trimmedSignal > min(window) & trimmedSignal < max(window)) | isnan(trimmedSignal)))));
	 catch Exception
		  beep
		  fprintf('%s\n',Exception.message);
		  trimTimeStampsError = MException('VerifyOutput:OutOfBounds',...
								'Results are outside the allowable limits');
		  throw(trimTimeStampsError);
	 end
end

function [spikesSignal] = test_spikeCounts(trimmedSignal,window,binWidth,expectedSpikes)
	spikesSignal = spikeCounts(trimmedSignal,window,binWidth);
	try
		  assert(all(all(spikesSignal == expectedSpikes)));
	 catch Exception
		  beep
		  fprintf('%s\n',Exception.message);
		  spikeCountsError = MException('VerifyOutput:OutOfBounds',...
								'Results are outside the allowable limits');
		  throw(spikeCountsError);
	 end
end

function [thisPsth] = test_psth(expectedSpikes,expectedPSTH)
	thisPsth = psth(expectedSpikes);
	try
		assert(all(thisPsth == expectedPSTH))
	catch Exception
		beep
		fprintf('%s\n',Exception.message);
		psthError = MException('VerifyOutput:OutOfBounds',...
						'Results are outside the allowable limits');
		throw(psthError);
	end
end

function test_significantSpan(vector, sig, thisSpan)
	try
		assert(all(significantSpan(vector, sig) == thisSpan))
	catch Exception
		beep
		fprintf('%s\n',Exception.message);
		sigSpanError = MException('VerifyOutput:OutOfBounds',...
						'Results are outside the allowable limits');
		throw(sigSpanError);
	end
end
