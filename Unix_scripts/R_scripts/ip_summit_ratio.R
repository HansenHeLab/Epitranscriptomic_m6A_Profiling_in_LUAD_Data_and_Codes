
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1];
ip_file <- args[2];
out_name <- args[3];

input <-read.table(input_file);

# for input 
reg_u <- sort(unique(input$V9));
cnt <- length(reg_u)
cov_input_list <- list()

for(i in 1:cnt) 
{
	idx <- reg_u[i] == input$V9;
	start_c <- input$V14[idx];
	end_c <- input$V15[idx];
	cov <- input$V16[idx];
	flag <- start_c[1]
	peak_left <- input$V2[idx][1];
	peak_right <- input$V3[idx][1];
	peak_l <- peak_right - peak_left;
	cov_v <- rep(0, peak_l);

	if(flag == -1)
	{
		cov_input_list[[i]] <- cov_v;
	}
	 else 
	{ 
		for (j in 1:length(cov))
		{
			start <- max(1, start_c[j]-peak_left+1);
			end <- min(peak_l, end_c[j]-peak_left);		# for bed file
			cov_v[start:end] =  as.numeric(as.character(cov[j]));
		}
		cov_input_list[[i]] <- cov_v;
	}
}

rm(input);

# for ip 
ip <- read.table(ip_file)
reg_u <- sort(unique(ip$V9));
cnt <- length(reg_u)
cov_ip_list <- list()

peak_left_v <-vector();
chr <- vector();
strand <- vector()

for(i in 1:cnt) 
{
	idx <- reg_u[i] == ip$V9;
	start_c <- ip$V14[idx];
	end_c <- ip$V15[idx];
	cov <- ip$V16[idx];
	flag <- start_c[1]
	peak_left <- ip$V2[idx][1];
	peak_right <- ip$V3[idx][1];
	peak_l <- peak_right - peak_left;
	cov_v <- rep(0, peak_l);
		
	chr[i] <- as.character(ip$V1[idx][1]);
	peak_left_v[i] <- peak_left;	
	strand[i] <- as.character(ip$V6[idx][1]);

	if(flag == -1)
	{
		cov_ip_list[[i]] <- cov_v;
	}
	 else 
	{ 
		for (j in 1:length(cov))
		{
			start <- max(1, start_c[j]-peak_left+1);
			end <- min(peak_l, end_c[j]-peak_left);		# for bed file
			cov_v[start:end] = as.numeric(as.character(cov[j]));
		}
		cov_ip_list[[i]] <- cov_v;
	}
}
rm(ip)

##

summit_pos = vector();
summit_fd = vector();

for (i in 1:cnt)
{
	fd_v <- cov_ip_list[[i]]/cov_input_list[[i]];
	fd_v[is.infinite(fd_v)] <- -1;			# #/0
	fd_v[is.nan(fd_v)] <- -2;  			# 0/0

	fd_max <- max(fd_v);
	idx_1 <- which(fd_v == fd_max);
     	ip_abs <- cov_ip_list[[i]][fd_v == fd_max]
        idx_2 <- ip_abs==max(ip_abs);
	idx <- idx_1[idx_2];

	idx_mid <- round((length(idx)+1)/2);
	offset <- idx[idx_mid];
        summit_fd[i] <- fd_max;
	summit_pos[i] <- peak_left_v[i]+offset-1;
} 

	

out <- cbind(chr,summit_pos, summit_pos+1, summit_fd,  reg_u, strand);
write.table(out, file=out_name, col.names=FALSE, row.name=FALSE, quote=FALSE, sep="\t");
