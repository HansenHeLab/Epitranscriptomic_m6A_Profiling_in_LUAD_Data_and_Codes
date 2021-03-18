
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1];
reg_len <- as.numeric(args[2]);
out_name <- args[3]

input <-read.table(input_file);

reg_u <- levels(input$V5);
cov_m <- matrix(0, length(input$V10), reg_len);

for(i in 1:length(reg_u)) 
{
	idx <- reg_u[i] == input$V5;
	cov <- input$V10[idx];
	strand <- input$V6[idx][1]
	start_c <- input$V8[idx];
	end_c <- input$V9[idx];
	peak_left <- input$V2[idx][1]

	for (j in 1:length(cov))                                                
    	 {                                                                       
             start <- max(1, start_c[j]-peak_left+1);                        
             end <- min(reg_len, end_c[j]-peak_left);         # for bed file  
             cov_m[i,start:end] = as.numeric(as.character(cov[j]));            
    	 } 

	if (strand == "-")
        {
		cov_m[i, ] <- rev(cov_m[i, ])
	} 
	                                                                      
		
}


idx_1 <- rowSums(cov_m) > 0;
cov_n <- cov_m[idx_1, ];
m_cov <- colMeans(cov_n);
colnames <- out_name;
write.table(m_cov, file=out_name, col.names=TRUE);


