function SPKobj = sz_spk_checkSPK(spk)
% Check the type of input.
% If string, try to load the file and check format,
% if matrix with spike times, use these to create a spike object,
% if struct ensure consistency of field names and return the struct.