#### cor plot function ###########################
# a : independent variable 
# b : dependent variable
# title:  title string
## fig_type :: png or pdf

cor.plot <- function(x, y, title, x_lab, y_lab, fig_type){
    
    ## png
    if (fig_type == "pdf" | fig_type == "PDF")
    {
        file_name <- paste(title, ".pdf", sep = "");
        pdf(file_name, width = 3.5, height = 3.5)
        par(mar = c(4, 3, 2, 0.5), mgp =c(2, 0.7, 0))
        #reg <- lm(y ~ x)
        cor <- cor.test(x , y)
        pval <- cor$p.value
        
        if(pval > 2.2e-16)
        {
            pval <- format(pval, scientific = T)
            sub_t <- paste("r = ", round(cor$estimate, 2), "; ", "p_value = ", pval, sep = "")        # correlation coefficient
            
        } else{
            
            sub_t <- paste("r = ", round(cor$estimate, 2), "; ", "p_value < 2.2e-16 ", sep = "") 
        } 
        
        plot(x, y, xlab = x_lab, ylab = y_lab, sub = title, main = sub_t, cex.main = 1, font.main = 3, cex.sub = 1.2, font.sub = 2 )
       #abline(b = reg$coefficients[2], a = reg$coefficients[1], lty = 2, col = "red" )
        abline(b = 1, a = 0, lty = 2, col = "gray" )
        dev.off()  
    } else {
    file_name <- paste(title, ".png", sep = "");
    png(file_name, width = 3.5, height = 3.5, units = "in", res = 300)
    par(mar = c(4, 3, 2, 0.5), mgp =c(2, 0.7, 0))
    #reg <- lm(y ~ x)
    cor <- cor.test(x , y)
    pval <- cor$p.value
    
    if(pval > 2.2e-16)
    {
        pval <- format(pval, scienctific = T)
        sub_t <- paste("r = ", round(cor$estimate, 2), "; ", "p_value = ", pval, sep = "")        # correlation coefficient
        
    } else{
        
        sub_t <- paste("r = ", round(cor$estimate, 2), "; ", "p_value < 2.2e-16 ", sep = "") 
    } 
        
    plot(x, y, xlab = x_lab, ylab = y_lab, sub = title, main = sub_t, cex.main = 1, font.main = 3, cex.sub = 1.2, font.sub = 2 )
    #abline(b = reg$coefficients[2], a = reg$coefficients[1], lty = 2, col = "red" )
    abline(b = 1, a = 0, lty = 2, col = "gray" )
    dev.off()
    }

}
