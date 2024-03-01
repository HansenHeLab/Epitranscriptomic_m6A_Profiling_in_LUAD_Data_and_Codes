create.subset.heatmap <- function(data) {
      
      key.min <- -4
      key.max <-  4
      key.colour.interval.num <- 50
      key.scale <- c(
            seq(key.min, 0, -key.min / key.colour.interval.num),
            seq(0, key.max, key.max / key.colour.interval.num)
      );
      key.scale <- unique(key.scale)
      
      rnaseq.heatmap <- create.heatmap(
            x = data,
            clustering.method = 'none',
            at = key.scale,
            colour.scheme = c('darkorchid4', 'white', 'darkgreen'),
            print.colour.key = FALSE,
            same.as.matrix = TRUE,
            scale.data = FALSE
      )
      
      return(rnaseq.heatmap)
}

create.covariate.heatmap <- function(annot) {
      numeric <- annot;

      numeric[which(annot$age <= 65),'age'] <- 'gray80'
      numeric[which(annot$age > 65),'age'] <- 'gray10'
      
      numeric['sex'] <- force.colour.scheme(annot$sex, scheme = "sex")
   
      numeric[which(annot$pathStage == "I"), 'pathStage'] <- '#FFF0E6'
      numeric[which(annot$pathStage == "II"), 'pathStage'] <- '#FFB380'
      numeric[which(annot$pathStage == "III" | annot$pathStage == "IV"), 'pathStage'] <- '#813D18'
      
      numeric[which(annot$smoking == "Ex-Smoker"), 'smoking'] <- '#5ab4ac'
      numeric[which(annot$smoking == "Never"), 'smoking'] <- '#f5f5f5'
      numeric[which(annot$smoking == "Current"), 'smoking'] <- '#d8b365'
      
      numeric[which(annot$mutation == "No"),'mutation'] <- 'white'
      numeric[which(annot$mutation == "No Data"),'mutation'] <- '#FF2323'
      numeric[which(annot$mutation == "KRAS"),'mutation'] <- '#7FC97F'
      numeric[which(annot$mutation == "EGFR"),'mutation'] <- '#386CB0'
      numeric[which(annot$mutation == "MET"),'mutation'] <- '#FFFF99'
      numeric[which(annot$mutation == "PIK3CA"),'mutation'] <- '#beaed4'
      
      heatmap.colours <- na.omit(unique(unlist(numeric)));
      for (i in 1:length(heatmap.colours)) {
            numeric[numeric == heatmap.colours[i]] <- i;
      }
      numeric <- sapply(numeric, as.numeric);
      
      cov.heatmap <- create.heatmap(
            x = numeric,
            clustering.method = 'none',
            same.as.matrix = TRUE,
            colour.scheme = heatmap.colours,
            grid.col = TRUE,
            grid.row = FALSE,
            total.colours = length(heatmap.colours) + 1,
            fill.colour = 'darkgray',
            print.colour.key = FALSE
      )
      
      return(cov.heatmap);
}

# set general variables
key.min <- -4
key.max <-  4
key.colour.interval.num <- 50;
key.scale <- c(
      seq(key.min, 0, -key.min / key.colour.interval.num),
      seq(0, key.max, key.max / key.colour.interval.num)
);
key.scale <- unique(key.scale);

covariates.legend <- list(
      legend = list(
            colours = c('gray80','gray10'),
            labels = c(expression(''<= 65),expression('' > 65)),
            title = expression(bold(underline('Age')))
      ),
      legend = list(
            colours = c('#FFB6C1', '#87CEFA'),
            labels = c('F', 'M'),
            title = expression(bold(underline('Sex')))
      ),
      legend = list(
            colours = c('#f5f5f5', '#5ab4ac','#d8b365'),
            labels = c('Never', 'Ex-Smoker','Current'),
            title = expression(bold(underline('Smoking')))
      ),
      legend = list(
            colours = c('#FFF0E6', '#FFB380','#813D18'), 
            labels = c('I', 'II','III and IV'),
            title = expression(bold(underline('Clinical Stage')))
      ),
      legend = list(
            colours = c('white', '#FF2323','#7FC97F','#386CB0','#FFFF99','#beaed4'),
            labels = c('No', 'No Data','KRAS','EGFR','MET','PIK3CA'),
            title = expression(bold(underline('Mutation')))
      ),
      
      legend = list(
            colours = c('darkorchid4', 'white', 'darkgreen'),
            labels = c(
                  as.character(round(key.scale[1], digits = 2)),
                  as.character(median(key.scale)),
                  as.character(round(key.scale[length(key.scale)], digits = 2))
            ),
            title = expression(bold(underline('Z'[Intensity]))),
            continuous = TRUE,
            tck = 0,
            at = c(0, 50, 100)
      )
)
