library(RColorBrewer)
library(pBrackets)

pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}

set.seed(6)

make.fig <- function(file.name) {

  g <- function() {

    layout(matrix(c(1,0,2,0,3,0,4,0,0,5,6,6,6,6,
                    1,0,2,0,3,0,4,0,0,5,0,0,0,0,
                    1,0,2,0,3,0,4,0,0,5,0,0,0,0,
                    1,0,2,0,3,0,4,0,0,5,0,0,0,0,
                    1,0,2,0,3,0,4,0,0,5,0,0,0,0),
                  5, 14, byrow=TRUE),
           heights=c(1,0.4,1.5,0.2,1.5),
           widths=c(0.1,0,0.25,0.05,0.25,0.05,0.25,0.05,0.05,0.25,0.3,0.07,0.07,0.07))
    par(oma=c(2,4,0,0), mar=c(0.1, 0.1, 1, 0.1),
        mgp=c(2,0.2,0), tcl=0, cex.axis=0.7, cex.main=0.7)
    cols <- c(gray(0.9), gray(0.7), brewer.pal(9, 'OrRd')[c(5,7)])

    visits <- 1:3
    n.visits <- length(visits)
    n.sites <- 20
    site.mat <- matrix(1,nrow=n.visits,ncol=n.sites)

    ## focal species range
    y.range <- 3:14
    n.sites.focal <- n.visits * length(y.range)

    make.panel <- function(yy) {
      ## detections of focal species
      focal.detections <- rbinom(n.sites.focal,size=1,prob=0.15)
      ## site visits
      visited <- rbinom(length(site.mat),size=1,prob=0.5)
      ## other species detected
      other.species <- rbinom(length(site.mat),size=1,prob=0.3)

      ## site.mat[visits,y.range] <- 2
      site.mat[visited==1] <- 2
      site.mat[other.species==1] <- 3
      site.mat[visits,y.range][focal.detections==1] <- 4
      ## site.mat[visits,y.site] <- rbinom(n.sites.focal,size=1,prob=0.2)+2

      image(site.mat, col=cols,
            xaxt='n', yaxt='n',
            font.main=1, main=sprintf('Era %d', yy))

      ## add black grid lines
      grid.x <- 1:(n.visits-1)/(n.visits-1) - 1/(2*(n.visits-1))

      grid.y <- c(min(y.range)-1,y.range)/(n.sites-1) - 1/(2*(n.sites-1))
      sapply(grid.x, function(x)
        lines(x=c(x,x),
              y=range(grid.y),
              col='yellow', lwd=0.75))
                                                      
      ## abline(v=grid.x, col=grey(0), lwd=0.75)
      abline(h=grid.y, col='yellow', lwd=0.75)
      site.mat
    }
    
    plot(NA, xlim=c(-0.005,0.005), ylim=c(0,1),
         xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
    brackets(0, 0.08, 0, 0.71, lwd=1)
    site.mat.1 <- make.panel(yy=1)
    site.mat.2 <- make.panel(yy=2)
    site.mat.2 <- make.panel(yy=3)
    mtext(expression('...'),
          side=1, line=-8.25, at=0.54, cex=1, adj=0, outer=TRUE)
    site.mat.2 <- make.panel(yy=10)
    
    plot(NA, xlim=c(0,1), ylim=c(0,1),
         xlab='', ylab='', xaxt='n', yaxt='n', bty='n')
    legend(x=0.04, y=1.1,
           legend=c('No visit occurred',
                    'Visit occurred, nothing detected',
                    'At least one species detected',
                    expression('Species'~italic(i)~'detected')),
           pch=rep(15, 4),
           col=cols,
           cex=0.7,
           pt.cex=1.7,
           bty='n')


    # ## SS_all,all
    # image(site.mat.1>0, xaxt='n', yaxt='n', col=c('black'))
    # image(site.mat.1>0, xaxt='n', yaxt='n', col=c('black'))
    # ## MS_range,all
    # site.mat.r <- site.mat.1*0
    # site.mat.r[,y.range] <- 1
    # image(site.mat.r, xaxt='n', yaxt='n', col=c('white','black'))
    # image(site.mat.r, xaxt='n', yaxt='n', col=c('white','black'))
    # ## MS_range,detected
    # site.mat.r.1 <- (site.mat.1==3 | site.mat.1==4)*site.mat.r
    # site.mat.r.2 <- (site.mat.2==3 | site.mat.2==4)*site.mat.r
    # image(site.mat.r.1, xaxt='n', yaxt='n', col=c('white','black'))
    # image(site.mat.r.2, xaxt='n', yaxt='n', col=c('white','black'))
    # ## MS_range,visits
    # site.mat.v.1 <- (site.mat.1==2 |
    #                  site.mat.1==3 |
    #                  site.mat.1==4)*site.mat.r
    # site.mat.v.2 <- (site.mat.2==2 |
    #                  site.mat.2==3 |
    #                  site.mat.2==4)*site.mat.r
    # image(site.mat.v.1, xaxt='n', yaxt='n', col=c('white','black'))
    # image(site.mat.v.2, xaxt='n', yaxt='n', col=c('white','black'))

    

    mtext(expression('Range of'),
          side=1, line=-8.25, at=-0.17, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('species'~italic(i)~','~italic(R[i])),
          side=1, line=-7.25, at=-0.17, cex=0.5, adj=0, outer=TRUE)
    
    mtext(expression('Site 1'),
          side=1, line=-18.2, at=-0.03, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('Site 2'),
          side=1, line=-17.2, at=-0.03, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('...'),
          side=1, line=-16.2, at=-0.03, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('Site 20'),
          side=1, line=-1.2, at=-0.04, cex=0.5, adj=0, outer=TRUE)

    v1 <- 0.07
    v3 <- 0.18
    v2 <- v1+(v3-v1)/2
    mtext(expression('1'),
          side=1, line=-0.1, at=v1, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('2'),
          side=1, line=-0.1, at=v2, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('3'),
          side=1, line=-0.1, at=v3, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('Interval'),
          side=1, line=0.8, at=0.1, cex=0.5, adj=0, outer=TRUE)

    shift <- 0.17
    mtext(expression('1'),
          side=1, line=-0.1, at=v1+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('2'),
          side=1, line=-0.1, at=v2+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('3'),
          side=1, line=-0.1, at=v3+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('Interval'),
          side=1, line=0.8, at=0.1+shift, cex=0.5, adj=0, outer=TRUE)
    
    shift <- 0.34
    mtext(expression('1'),
          side=1, line=-0.1, at=v1+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('2'),
          side=1, line=-0.1, at=v2+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('3'),
          side=1, line=-0.1, at=v3+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('Interval'),
          side=1, line=0.8, at=0.1+shift, cex=0.5, adj=0, outer=TRUE)
    
    shift <- 0.51
    mtext(expression('1'),
          side=1, line=-0.1, at=v1+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('2'),
          side=1, line=-0.1, at=v2+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('3'),
          side=1, line=-0.1, at=v3+shift, cex=0.5, adj=0, outer=TRUE)
    mtext(expression('Interval'),
          side=1, line=0.8, at=0.1+shift, cex=0.5, adj=0, outer=TRUE)

    # mtext(expression(~italic(SS['all,all'])~', '~italic(MS['all,all'])),
    #       side=1, line=-13.3, at=0.44, cex=0.5, adj=0, outer=TRUE)
    # mtext(expression(~italic(MS['range,all'])),
    #       side=1, line=-13.3, at=0.71, cex=0.5, adj=0, outer=TRUE)
    # mtext(expression(~italic(MS['range,detected'])),
    #       side=1, line=-6.3, at=0.47, cex=0.5, adj=0, outer=TRUE)
    # mtext(expression(~italic(MS['range,visits'])),
    #       side=1, line=-6.3, at=0.71, cex=0.5, adj=0, outer=TRUE)
    
  }
  
  pdf.f(g, file = file.name, height=2.75, width=4.5)
}

make.fig('figures/model_v4.pdf')

