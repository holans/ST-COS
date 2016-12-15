function fid = openlogger(filename)
	global log_fid
	log_fid = fopen(filename,'w');
	fid = log_fid;
end

