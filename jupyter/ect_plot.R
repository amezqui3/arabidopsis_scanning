library(ggplot2)
library(colorRamps)
library(viridis)
library(ggthemes)

setwd('/home/ejam/documents/arabidopsis_scanning/petiole/ecc')

sname <- 'col-0_pot0_leaf0'
axis <- '2'
TT <- 32

verts <- as.matrix(read.csv(paste(sname,'_coords.csv',sep = ''), header=F))
colrs <- unlist(as.vector(read.csv(paste(sname,'_vals_',axis,'.csv',sep = ''), header=F)))
sigma <- 1.05*unlist(as.vector(read.csv(paste(sname,'_sigma.csv',sep = ''), header=F)))
col <- (viridis::plasma(TT)) 

lout <- 50
MX <- cbind(seq(-sigma[1], sigma[1], length.out = lout), rep.int(0,lout), rep.int(0,lout))
MY <- cbind(rep.int(0,lout), seq(-sigma[2], sigma[2], length.out = lout), rep.int(0,lout))
MZ <- cbind(rep.int(0,lout), rep.int(0,lout), seq(-sigma[3], sigma[3], length.out = lout))

rgl::open3d()
rgl::par3d(windowRect = c(50,50,1000,455))
bar <- rbind(apply(verts, 2, max), apply(verts, 2, min))
rgl::plot3d(bar, xlab = ' ', ylab = ' ', zlab = ' ',
            col = 'white', type = 'n', lwd = 4, radius = rd, size = 0.1, alpha=0,
            ann=FALSE, axes=FALSE)
rgl::rgl.viewpoint(userMatrix = um2,zoom = 0.31)
rgl::aspect3d('iso')

#rgl::points3d(verts, col=col[colrs], sz=1)
{
  rgl::segments3d(MX,
                  lwd=10, col = '#dc3220')
  rgl::segments3d(MY,
                  lwd=8, col = '#005ab5')
  rgl::segments3d(c(0,0),c(0,0), c(sigma[3],-sigma[3]),
                  lwd=6, col = '#00e00c')
}
#m2 <- rgl::par3d()$userMatrix
sz <- 3
i<- 0
rgl::rgl.snapshot(paste('ecc_axis_',axis,'_',formatC(i, digits=1, format='d', flag='0'),'.png', sep=''), fmt='png')
for( i in 1:length(col)){
  rgl::rgl.viewpoint(userMatrix = um2,zoom = 0.31)
  rgl::points3d(verts[which(colrs <= i),],col =col[colrs[which(colrs <= i)]], size = sz)
  rgl::rgl.snapshot(paste('ecc_axis_',axis,'_',formatC(i, digits=1, format='d', flag='0'),'.png', sep=''), fmt='png')
}

#################################################################################################################3
#################################################################################################################3
#################################################################################################################3

setwd('/home/ejam/documents/arabidopsis_scanning/petiole/ecc')

axis <- '2'
ecc <- ceiling(0.51*as.data.frame(read.csv(paste(sname,'_ecc_',axis,'.csv',sep = ''), header=F)))
ecc$Index <- 1:dim(ecc)[1]
col <- (viridis::magma(dim(ecc)[1])) 
i <- 32
p <- ggplot(ecc, aes(x=Index, y=V1))+
  geom_path()+
  geom_point(shape=21, col='black', size=3, fill=col[1:i], stroke=1) +
  geom_point(aes(x=1, y=min(V1)), colour="blue", alpha = 0) +
  geom_point(aes(x=length(V1), y=max(V1)), colour="blue", alpha = 0) +
  #scale_x_discrete(labels=c('2'='2','4'='4','6'='6','14'='14','18'='18','32'='32')) +
  scale_x_continuous(labels=c(expression('t'[2]),expression('t'[6]),expression('t'[8]),
                              expression('t'[12]),expression('t'[20]),expression('t'[32])), 
                     breaks=c(2,6,8,12,20,32)) +
  theme_classic() +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title = element_text(size=16),
        axis.text = element_text(size=16),
        axis.title.y = element_text(angle=0, vjust=0.5),
        panel.grid.major.x = element_line(color = 'grey50', linetype=2)) +
  labs(x = expression('threshold'),
       y = expression(paste(chi)),
       title = paste('ECC (Z-axis direction, T =',TT,')')) +
  geom_rangeframe() 
p
filename <- paste(sname,'_ecc_Z_',formatC(i, digits = 1, format='d', flag='0'), sep='')
ggsave(paste(filename,'.png', sep=''), plot=p, device='png',
       width = 141.1, height = 105.7, units = 'mm', dpi=72)
ggsave(paste(filename,'.pdf', sep=''), plot=p, device='pdf',
       width = 141.1, height = 105.7, units = 'mm', dpi=72)

#################################################################################################################3
#################################################################################################################3
#################################################################################################################3

setwd('/home/ejam/documents/arabidopsis_scanning/petiole/ect')

ecc <- as.data.frame(read.csv(paste(sname,'_ecc_',axis,'.csv',sep = ''), header=F))
ecc$Index <- 1:dim(ecc)[1]
col <- (viridis::magma(dim(ecc)[1])) 
ect <- ceiling(0.51*as.data.frame(read.csv(paste(sname,'_ect.csv',sep = ''), header=F)))
ect$Index <- 1:dim(ect)[1]


df <- data.frame( 
  X = 1:length(ect), 
  ECT = ect
)

i <- 90
cols <- rep(col,3)

title <- paste('Euler Characteristic Transform (3 directions,',TT,'thresholds each)')

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
  return(metadir)
}

for(i in 1:nrow(ect)){
  metadir <- def_metadir(i, TT)
  p <- ggplot(ect[1:i,]) +
    geom_rect(aes(xmin=start, xmax=end, fill=Direction),
              ymin=-Inf, ymax=Inf, alpha=0.3, data=metadir) +
    geom_path(aes(Index, V1))+
    geom_point(aes(Index, V1),shape=21, size=3, fill=cols[1:i], stroke=1) +
    geom_point(aes(x=1, y=min(ect)), colour="blue", alpha = 0) +
    geom_point(aes(x=nrow(ect), y=max(ect[,1])), colour="blue", alpha = 0) +
    scale_x_continuous(breaks=8*(0:12)) +
    theme(plot.title = element_text(size=30, hjust = 0.5),
          axis.title.y = element_text(size=30, angle=0, vjust=0.5),
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
    scale_fill_manual(values=c('#dc3220','#005ab5','#00e00c'))
  p
  filename <- paste(sname,'_ect_',formatC(i, digits = 1, format='d', flag='0'), sep='')
  w <- 3*146
  h <- 1.2*93
  ggsave(paste(filename, '.pdf', sep=''), plot=p, device='pdf', width = w, height = h, units='mm')
  ggsave(paste(filename, '.jpg', sep=''), plot=p, device='jpg', width = w, height = h, units = 'mm', dpi=96)
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
ggsave(paste(filename, '.pdf', sep=''), plot=p, device='pdf', width = w, height = h, units='mm')
ggsave(paste(filename, '.jpg', sep=''), plot=p, device='jpg', width = w, height = h, units = 'mm', dpi=96)
