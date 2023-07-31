.pairsCorColor <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {

    # set P value cutoff
    pvalCo = 0.05

    # get the x and y data to use the other code
    x <- GGally::eval_data_col(data, mapping$x)
    y <- GGally::eval_data_col(data, mapping$y)

    # calc pairwise pearson correlation
    ct <- corrplot::cor.test(x,y, method = "pearson")

    # get P value
    pval = ct$p.value

    # get correlation coefficient
    r <- unname(ct$estimate)
    rt <- format(r, digits=3)[1]
    tt <- as.character(rt)

    # plot the cor value
    p <- GGally::ggally_text(
        label = tt,
        mapping = aes(),
        xP = 0.5, yP = 0.5,
        size = 4,
        color = color,
        ...
    ) +
        theme_void() +
        theme(
            panel.background=element_rect(fill="white"),
            panel.grid.minor=element_blank(),
            panel.grid.major=element_blank()
        )

    # set color scheme for background
    colfunc = colorRampPalette(c("#cfd6e2", "#303c50"))
    corColors = colfunc(n = 10)

    # apply color scheme for background
    if(r >= .9){
        corCol = corColors[10]
    }
    else if(r >= .8){
        corCol = corColors[9]
    } else if(r >= .7){
        corCol = corColors[8]
    } else if(r >= .6){
        corCol = corColors[7]
    } else if(r >= .5){
        corCol = corColors[6]
    } else if(r >= .4){
        corCol = corColors[5]
    } else if(r >= .3){
        corCol = corColors[4]
    } else if(r >= .2){
        corCol = corColors[3]
    } else if(r >= .1){
        corCol = corColors[2]
    }
    else {
        corCol = "white"
    }

    # set background color to grey for not significant correlations
    if (pval >= pvalCo) {
        corCol = "#A9A9A9"
    }

    p <- p + theme(
        panel.background = element_rect(fill = corCol,
                                        color = corCol)
    )
    return(p)
}

.pairsDensity <- function(data, mapping, ...){
    # init local variables
    y <- NULL

    data = data %>% drop_na()
    # compute quantiles to display
    x <- GGally::eval_data_col(data, mapping$x)
    x = as.numeric(x)
    dt = data.frame(x = c(1:length(x)), y = x)
    dens <- density(dt$y)
    df <- data.frame(x=dens$x, y=dens$y)
    probs <- c(0, 0.5, 1)
    quantiles <- quantile(dt$y, prob=probs)
    quant <- factor(findInterval(df$x,quantiles))
    df$quant <- factor(findInterval(df$x,quantiles))

    # make plot
    p = ggplot(df, aes(x,y)) +
        geom_line() + geom_ribbon(aes(ymin=0, ymax=y), fill = "#ECE5C7") +
        scale_x_continuous(breaks=quantiles) +
        scale_fill_brewer(guide="none") +
        theme_bw() +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
    return(p)
}

.parisCorrelationPlot <- function(data, mapping, max.value) {

    X <- GGally::eval_data_col(data, mapping$x)
    Y <- GGally::eval_data_col(data, mapping$y)

    p = ggplot(data, mapping) +
        geom_bin2d(bins = 100) +
        viridis::scale_fill_viridis() +
        theme_bw() +
        theme(legend.position = "none") +
        scale_x_continuous(limits = c(-1, ceiling(max.value))) +
        scale_y_continuous(limits = c(-1, ceiling(max.value)))

    return(p)
}
