plotnetwork <-
function(datainput,                    # The correlation matrix, usually a lower matrix
         n_levels = 6,                 # Number of types of segments to show the correlation, note interval should <12
         xlim = c(-2.5, 5),            # horizontal range of the canvas
         ylim = c(-2.5, 2.5),          # vertical range of the canvas
         node_size = 3,                # size of label
         lwd.var = TRUE,               # logical, if the segments width should vary with the absolute value
         lwd = 4,                      # width of the segments for connecting the sites, default 1
         label_dist = 1.2,             # Distance of text labels from each corner.
         show.node = TRUE,             # Whether the nodes in the figure should be labeled
         show.text.label = TRUE,       # Whether the text label should be drawn.
         linecol = c("orange", "blue"),# Colours showing positve or negative correlation.
         # The lines representing positive correlations use the first element
         show.legend = TRUE,           # If the legend should be shown
         valuename = "r",              # Name of the variable shown in the legend
         legendx = 3,                  # the starting position of the legend x
         legendy = -2,                 # the starting position of the legend y
         legend_line_space = 1,        # the space between the lines in the legend
         legend_linelength = 0.3,      # length of the segment in the legend
         adjust_legend_x = 0,          # adjusting the position of legend (x)
         adjust_legend_y = 0,          # adjusting the position of legend (y)
         digits = 2,                   # Number of digits shown in the legend, 
         # it will be used in generating the breaks as well
         ...                           # other parameters related with plot.default
        ){

    datainput <- as.matrix(datainput)

    if(any(is.na(datainput))){
        stop("NA is not allowed")
    }
    if(n_levels > 6|n_levels < 2){
        stop("Number of levels should be between 2 and 6")
    }
    if (!nrow(datainput)==ncol(datainput)){
        stop("The input is not a correlative matrix.")
    }
    npoints = nrow(datainput)
    if (npoints > 10){
        warning("Too many sites to show! Please use matrix with fewer (e.g., <10) dimensions")
    }

    plot(0,0, xlim = xlim, ylim = ylim,
         type = "n", axes = FALSE, xlab="", ylab="", ...)

    seqn  <- (seq(0, n_levels^2, by = n_levels))/(n_levels^2)
    limit0 <- sort(c(min(datainput) + (max(datainput)-min(datainput))*seqn))

    subangle <- 2*pi/npoints

    r=2

    xdep <- c()
    ydep <- c()
    txdep <- c()
    tydep <- c()

    for (k in 1:npoints){
        r = r
        x = r*sin(subangle * k)
        y = r*cos(subangle * k)
        tx = label_dist*x
        ty = label_dist*y
        xdep <- append(xdep, x)
        ydep <- append(ydep, y)
        txdep <- append(txdep, tx)   # Position for labels x
        tydep <- append(tydep, ty)   # Position for labels y
    }

    ltys <- sort(1:n_levels, decreasing = TRUE)

    for (i in 1:length(xdep)){
        for (j in 1:i){ # Only loop through the lower matrix
            for (m in 1:n_levels){
                colori <- ifelse(datainput[i,j] > 0, linecol[1], linecol[2])
                    if (datainput[i,j] >= limit0[m] & datainput[i,j] < limit0[m + 1]){
                        segments(x0 = xdep[i],
                                 y0 = ydep[i],
                                 x1 = xdep[j],
                                 y1 = ydep[j],
                                 lty = ltys[m],  # m is the variable telling the what kind of lines to draw
                                 col= colori,
                                 lwd = ifelse(lwd.var, abs(mean(limit0[m:(m+1)]))*lwd, lwd))
                 } else {
                 if(datainput[i,j] == limit0[m + 1]){
                         segments(x0 = xdep[i],
                                 y0 = ydep[i],
                                 x1 = xdep[j],
                                 y1 = ydep[j],
                                 lty = ltys[m],  # m is the variable telling the what kind of lines to draw
                                 col= colori,
                                 lwd = ifelse(lwd.var, abs(mean(limit0[m:(m+1)]))*lwd, lwd))
                    }
                }
            }
        }
    }

   if(show.node){
       points(xdep, ydep, cex=node_size, pch = 21, bg = "white", lwd = lwd/2)
       text(xdep, ydep, 1:npoints, cex = node_size/4)
   }
   if(show.text.label){
       text(txdep, tydep, colnames(datainput))
   }
   if(show.legend){
       for (n in 1:n_levels){
            col_n <- ifelse(mean(limit0[n] + limit0[n+1]) > 0, linecol[1], linecol[2])

            segments(legendx + adjust_legend_x,
                     legendy + legend_line_space * n * 0.25 + adjust_legend_y,
                     legendx + adjust_legend_x + legend_linelength,
                     legendy + legend_line_space * n * 0.25 + adjust_legend_y,
                     lty = ltys[n],
                     #col = col_n,
                     lwd = ifelse(lwd.var,
                                  ifelse(abs(mean(limit0[n:(n+1)]))*lwd < 0.01,
                                              0.1,
                                              abs(mean(limit0[n:(n+1)]))*lwd),
                                              lwd)
                     )
            if(n < length(limit0)){

                    text(legendx + adjust_legend_x + 5/4 + 0.5 + legend_linelength, 
                     legendy + legend_line_space*n*0.25 + adjust_legend_y, 
                     paste(formatC(sprintf("%.2f", limit0[n+1]), width = 5)))

                    text(legendx + adjust_legend_x + 4/4 + 0.5 + legend_linelength, 
                     legendy + legend_line_space*n*0.25 + adjust_legend_y, 
                     expression(""<""))

                    text(legendx + adjust_legend_x + 3/4 + 0.5 + legend_linelength, 
                     legendy + legend_line_space*n*0.25 + adjust_legend_y, 
                     valuename)


                    text(legendx + adjust_legend_x + 2/4 + 0.5 + legend_linelength, 
                    legendy + legend_line_space*n*0.25 + adjust_legend_y, 
                    expression(""<=""))

                    text(legendx + adjust_legend_x + 1/4 + 0.5 + legend_linelength, 
                    legendy + legend_line_space*n*0.25 + adjust_legend_y, 
                    paste(formatC(sprintf("%.2f", limit0[n]), width = 5)))

            }
          }
    }
}
