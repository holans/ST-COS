function fid = openlogger(filename)
	global log_fid
	log_fid = fopen(filename,'a');
	fid = log_fid;
end

