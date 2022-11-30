#!/usr/bin/env Rscript
# this script does what?
library(Polychrome)
library(ggplot2)
library(Cairo)

# Examples of meth beta values
# alternating intensities Bx (odd idx) BrM (even)
c1 <- c(0.5093399, 0.5913508, 0.5904641, 0.2042697, 0.8306049, 0.6609125, 0.2571777, 0.2317477, 0.2320412, 0.5177484, 0.5879986, 0.3454807, 0.4236549, 0.2134360, 0.2001799, 0.1690186, 0.4541413, 0.4270391, 0.2607895, 0.2271527, 0.3270316, 0.1473630, 0.6041846, 0.4239043, 0.3024745, 0.2357595, 0.3710093, 0.1859094, 0.6320090, 0.4505895, 0.4469454, 0.3601695, 0.4710973, 0.2300328)
# means was "0.4412437", "0.3306991")
                                                                                                                                   
c2 <- c(0.37519869, 0.36073997, 0.63703805, 0.29571714, 0.79699470, 0.56506777, 0.36893794, 0.27713252, 0.16874741, 0.41039883, 0.59069011, 0.25808665, 0.39241958, 0.11808662, 0.20107892, 0.22469188, 0.40102529, 0.46717153, 0.23636329, 0.23328335, 0.18070411, 0.06673689, 0.54374672, 0.11598903, 0.23863292, 0.20189831, 0.50014710, 0.51181998, 0.53060862, 0.44849552, 0.64114855, 0.54272941, 0.42291618, 0.16834783)
# means "0.42508225", "0.30978784")
                                                                                                                                   
c3 <- c(0.6308576, 0.7854971, 0.5044553, 0.6522234, 0.3659531, 0.3672468, 0.3604295, 0.7106267, 0.5513593, 0.5806608, 0.5243173, 0.6659574, 0.4585073, 0.5550123, 0.5875399, 0.6665482, 0.5138795, 0.4956739, 0.5448193, 0.5383910, 0.3628784, 0.4784778, 0.4328459, 0.5613406, 0.6516721, 0.5529366, 0.4590554, 0.3905431, 0.5629083, 0.6837479, 0.5676619, 0.4738731, 0.6158811, 0.5860428)
# means "0.5114718", "0.5732235")

c1 <- c1[7:12]
c2 <- c2[7:12]
c3 <- c3[7:12]

tlab <-c("first", "sec", "third")

ll <- length(c1)
c1d <- c1[seq(2,ll,2)] - c1[seq(1,ll,2)]
c2d <- c2[seq(2,ll,2)] - c2[seq(1,ll,2)]
c3d <- c3[seq(2,ll,2)] - c3[seq(1,ll,2)]
ll2 <- as.integer(ll/2)
pal2 <- alphabet.colors(ll2)
names(pal2) <- NULL # this is ultra important ... scale_fill_manual's major quirk is gettign confused by names stringvectors
# quickly set alpha
pal2 <- paste0(pal2, "99")
mypal <- c("#25739155", "#dd22aa55", "#11ee6655")

df <- data.frame(Diffs=round(c(c1d,c2d,c3d), digits=3), Cg=factor(rep(paste0("cg", 1:3), each=ll2)), Pairs=factor(rep(paste0("P", 1:ll2), 3)))

CairoPNG("ctri.png", 800, 800)
gg <- ggplot(df, aes(x=Cg, y=Diffs, fill=Pairs, label=Diffs)) + 
      # geom_bar(stat="identity", aes(fill=pal2)) +
      # geom_bar(stat="identity") +
      geom_col() +
      scale_fill_manual(values=pal2) +
      geom_text(size=3, position = position_stack(vjust = 0.5)) +
      geom_text(tlab, aes(label=tlab), position=position_dodge(width=0.9), vjust=-0.25) +
      # scale_fill_manual(values=c("#25739155", "#dd22aa55", "#11ee6655")) +
      geom_hline(yintercept=0) +
      theme_light() +
      theme(axis.line.y=element_line())
show(gg)
dev.off()
