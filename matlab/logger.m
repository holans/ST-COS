function logger(msg)
	global log_fid
	fprintf(log_fid, '%s - ', datestr(now, 31));
	fprintf(log_fid, '%s', msg);
	fprintf(log_fid, '\n');
end

