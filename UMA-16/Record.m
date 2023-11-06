fs = 16000; %Setting the sampling frequency
audioFrameLength = 1024; %Setting the buffer size
%Calling the UMA card
deviceReader = audioDeviceReader(...
 'Device', 'miniDSP ASIO Driver',...
 'Driver', 'ASIO', ...
 'SampleRate', fs, ...
 'NumChannels', 16,...
 'OutputDataType','double',...
 'SamplesPerFrame', audioFrameLength);
setup(deviceReader)
fileWriter = dsp.AudioFileWriter('PFE Test.wav','FileFormat','WAV', 'DataType', 'int16');
disp('Speak into microphone now.')
%Starting recording on tic
tic
while toc< 20
    acquiredAudio = step(deviceReader);
    step(fileWriter, acquiredAudio);
end
release(deviceReader);
release(fileWriter);
%End of recording
disp('Recording complete.')

