
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1];
ip_file <- args[2];
out_name <- args[3];
gf_file <- args[4];

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
gene <- vector();
strand <- vector();

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
		
	gene[i] <- as.character(ip$V4[idx][1]);
	strand[i] <- as.character(ip$V6[idx][1]);
	peak_left_v[i] <- peak_left;	

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

peak_ip_pileup <- vector()
peak_input_pileup <- vector()
peak_fd = vector();

for (i in 1:cnt)
{
	cov_ip_m <- mean(cov_ip_list[[i]])
        cov_input_m <- mean(cov_input_list[[i]])
        
	peak_ip_pileup[i] <- cov_ip_m;
	peak_input_pileup[i] <- cov_input_m;
	peak_fd[i] <- cov_ip_m/cov_input_m	
} 

	

out <- as.data.frame(cbind(gene,peak_input_pileup, peak_ip_pileup, peak_fd,  reg_u, strand));

gf_in <- paste("cut -f5,7 ", gf_file)   # read in genomefeatures infomation
gf <- read.table(pipe(gf_in))
for (k in 1: length(out$reg_u))
	{
		idx_2 <- out$reg_u[k]==gf$V1
		out$V7[k] <- as.character(gf$V2[idx_2][1])                                                                                            
	}


write.table(out, file=out_name, col.names=FALSE, row.name=FALSE, quote=FALSE, sep="\t");
