library(rgl) # needed for 3d graphics
library(ggplot2)
setwd('/home/ejam/documents/arabidopsis_scanning/pots')

bname <- 'col-0'
filename <- 'col-0_pot1_incylinder.csv'
verts <- as.matrix(read.csv(paste(bname,'/',filename, sep=''), header=FALSE))
verts <- verts[, c(3,1,2)]
sverts <- verts[1:(nrow(verts)%/%2), ]

rgl::open3d()
rgl::par3d(windowRect = c(50,50,585,455))
bar <- rbind(apply(verts, 2, max), apply(verts, 2, min))
rgl::plot3d(bar, xlab = ' ', ylab = ' ', zlab = ' ',
            col = 'white', type = 'n', lwd = 4, radius = rd, size = 0.1, alpha=0,
            ann=FALSE, axes=FALSE)
rgl::aspect3d('iso')

rgl::points3d(verts, col='forestgreen', sz=0.1)


degs <- seq(0, 360, by=10)
degs2 <- seq(-10, 10, length.out = length(degs)%/%2 + 1)
degs3 <- seq(10, -10, length.out = length(degs2))
degs4 <- c(degs2, degs3[-1])

i <- 1
for(i in 1:length(degs4)){
  rgl::rgl.viewpoint(degs[i],degs4[i],zoom = 0.28)
  rgl::rgl.snapshot(paste('col0_pot1_',formatC(i, digits=1, format='d', flag='0'),'.png', sep=''), fmt='png')
}
  


rgl::plot3d(verts , col=col[colrs], size = 12, 
            xlab='',ylab='',zlab='',
            ann=FALSE, axes=FALSE)
rgl::rgl.viewpoint(10,-10,zoom = 0.3)
rgl::aspect3d('iso')
rgl::rgl.snapshot(paste('height/height_filtered.png', sep=''), fmt='png')

#################################################################################################################3
#################################################################################################################3
#################################################################################################################3

rgl::open3d()
rgl::par3d(windowRect = c(50,50,900,700))
bar <- rbind(apply(verts, 2, max), apply(verts, 2, min))
rgl::plot3d(bar, xlab = '', ylab = '', zlab = '',
            col = 'white', type = 'n', lwd = 4, radius = rd, size = 0.1, alpha=0,
            axes=FALSE, ann=FALSE)
rgl::rgl.viewpoint(150,27,zoom = 0.58)
rgl::aspect3d('iso')
sz <- 1
rgl::points3d(verts, col='grey', size=sz)
i <- 0

for( i in 1:length(col)){
  sz <- 8 + sz
  rgl::rgl.viewpoint(150,27,zoom = 0.58)
  rgl::points3d(verts, col=col[i], size = sz)
  rgl::rgl.snapshot(paste('vrips/vrips_ver8_filter_',formatC(i, digits=1, format='d', flag='0'),'.png', sep=''), fmt='png')
}

rgl::rgl.viewpoint(155,25,zoom = 0.46)

#################################################################################################################3
#################################################################################################################3
#################################################################################################################3

sname <- 'seed_10_0'
axis <- 'Z'
TT <- 32

ecc <- ceiling(0.51*as.data.frame(read.csv(paste(sname,'_coords_ecc_',axis,'.csv',sep = ''), header=F)))
ecc$Index <- 1:dim(ecc)[1]


ect <- ceiling(0.51*as.data.frame(read.csv(paste(sname,'_coords_ect_T',TT,'.csv',sep = ''), header=F)))
ect$Index <- 1:dim(ect)[1]

i <- 20
col <- (viridis::magma(TT))
cols <- rep(col,3)

dst <- '~/documents/barley_stacks/manuscript2/ecc/'
title <- paste('ECC (',axis,'-axis direction , ',TT,' thresholds)', sep='')

for(i in 1:nrow(ecc)){
  metadir <- def_metadir(i, TT) 
  p <- ggplot(ecc[1:i,]) +
    geom_rect(aes(xmin=start, xmax=end), fill='green3',
              ymin=-Inf, ymax=Inf, alpha=0.15, data=metadir) +
    geom_path(aes(Index, V1))+
    geom_point(aes(Index, V1), shape=21, col='black', size=3, fill=cols[1:i], stroke=1) +
    geom_point(aes(x=1, y=min(ecc)), colour="blue", alpha = 0) +
    geom_point(aes(x=nrow(ecc), y=max(ecc[,1])), colour="blue", alpha = 0) +
    #scale_x_discrete(labels=c('2'='2','4'='4','6'='6','14'='14','18'='18','32'='32')) +
    #scale_x_continuous(labels=c(expression('t'[2]),expression('t'[4]),expression('t'[6]),
    #                            expression('t'[8]),expression('t'[14]),expression('t'[32])), 
    #                   breaks=c(2,4,6,8,14,32)) +
    #scale_x_continuous(breaks=4*(0:8)) +
    theme(plot.title = element_text(size=25, hjust = 0.5),
          axis.title.y = element_text(size=1.5*20),
          axis.title.x = element_text(size=20),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          #panel.grid.major.x = element_line(color = 'red', linetype=2)
    ) +
    labs(x = 'threshold',
         y = expression(paste(chi)),
         title = title)
  #+ geom_rangeframe() 
  #p
  filename <- paste(sname,'_ecc_',formatC(i, digits = 1, format='d', flag='0'), sep='')
  #ggsave(paste(dst,filename, '.pdf', sep=''), plot=p, device='pdf', width = 1*146, height = 1*93, units='mm')
  ggsave(paste(dst,filename, '.jpg', sep=''), plot=p, device='jpg', width = 1.25*146, height = 93, units = 'mm', dpi=96)
}

i <- 0
p <- ggplot() +
  geom_point(aes(x=1, y=min(ecc)), colour="blue", alpha = 0) +
  geom_point(aes(x=nrow(ecc), y=max(ecc[,1])), colour="blue", alpha = 0) +
  theme(plot.title = element_text(size=25, hjust = 0.5),
        axis.title.y = element_text(size=1.5*20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
  ) +
  labs(x = 'threshold',
       y = expression(paste(chi)),
       title = title)
p
filename <- paste(sname,'_ecc_',formatC(i, digits = 1, format='d', flag='0'), sep='')
#ggsave(paste(dst,filename, '.pdf', sep=''), plot=p, device='pdf', width = 1*146, height = 1*93, units='mm')
ggsave(paste(dst,filename, '.jpg', sep=''), plot=p, device='jpg', width = 1.25*146, height = 93, units = 'mm', dpi=96)
p

filename <- paste(sname,'_ect_',formatC(i, digits = 1, format='d', flag='0'), sep='')
ggsave(paste(filename, '.pdf', sep=''), plot=p, device='pdf', width = 3*146, height = 1*93, units='mm')
ggsave(paste(filename, '.pdf', sep=''), plot=p, device='pdf', width = 1.5*146, height = 2*93, units='mm')

#################################################################################################################3
#################################################################################################################3
#################################################################################################################3

title <- paste('Euler Characteristic Transform (3 directions,',TT,'thresholds each)')
i <- 24

def_metadir <- function(i, TT){
  start <- c(TT*(0:2)+1)
  end <- TT*(1:3)
  Direction = c('X','Y','Z')
  
  top <- ceiling(i/TT)
  xstart <- start[1:top]
  xend <- end[1:top]
  xend[top] <- i
  metadir <- data.frame(
    start = xstart,
    end = xend,
    Direction = Direction[1:top]
  )
}

metadir <- def_metadir(i, TT)

i <- 49

for(i in 1:nrow(ect)){
  metadir <- def_metadir(i, TT)
  p <- ggplot(ect[1:i,]) +
    geom_rect(aes(xmin=start, xmax=end, fill=Direction),
              ymin=-Inf, ymax=Inf, alpha=0.15, data=metadir) +
    geom_path(aes(Index, V1))+
    geom_point(aes(Index, V1),shape=21, size=3, fill=cols[1:i], stroke=1) +
    geom_point(aes(x=1, y=min(ect)), colour="blue", alpha = 0) +
    geom_point(aes(x=nrow(ect), y=max(ect[,1])), colour="blue", alpha = 0) +
    scale_x_continuous(breaks=8*(0:12)) +
    theme(plot.title = element_text(size=30, hjust = 0.5),
          axis.title.y = element_text(size=30),
          axis.title.x = element_text(size=20),
          axis.text = element_text(size=20),
          legend.position = 'bottom',
          legend.text = element_text(size=20),
          legend.title = element_text(size=20),
          legend.background = element_rect(color='grey', fill='grey95')
    ) +
    labs(x = 'threshold',
         y = expression(paste(chi)),
         title = title) +
    scale_fill_manual(values=c('red','blue','green3'))
  p
  filename <- paste(sname,'_ect_',formatC(i, digits = 1, format='d', flag='0'), sep='')
  ggsave(paste(dst, filename, '.pdf', sep=''), plot=p, device='pdf', width = 3*146, height = 1.2*93, units='mm')
  #ggsave(paste(filename, '.jpg', sep=''), plot=p, device='jpg', width = 146, height = 93, units = 'mm', dpi=96)
}

i<-0
p <- ggplot()+
  geom_point(aes(x=c(1,nrow(ect)), y=c(min(ect),max(ect[,1]))), colour="blue", alpha = 0) +
  scale_x_continuous(breaks=8*(0:12)) +
  theme(plot.title = element_text(size=30, hjust = 0.5),
        axis.title.y = element_text(size=30),
        axis.title.x = element_text(size=20),
        axis.text = element_text(size=20)
  ) +
  labs(x = 'threshold',
       y = expression(paste(chi)),
       title = title)
p
filename <- paste(sname,'_ect_',formatC(i, digits = 1, format='d', flag='0'), sep='')
ggsave(paste(dst, filename, '.pdf', sep=''), plot=p, device='pdf', width = 3*146, height = 1.2*93, units='mm')
#ggsave(paste(filename, '.jpg', sep=''), plot=p, device='jpg', width = 146, height = 93, units = 'mm', dpi=96)