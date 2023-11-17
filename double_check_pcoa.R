ggplot(data = subset(mcav.pca, habitat %in% "nearshore")) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  #scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
  #                  values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))

ggplot(data = subset(mcav.pca, habitat %in% "offshore")) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  #scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
  #                  values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))
ggplot(data = subset(mcav.pca, habitat %in% "deep")) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  #scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
   #                 values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))



ggplot(data = subset(ssid.pca, habitat %in% "nearshore")) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  #scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
  #                  values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))

ggplot(data = subset(ssid.pca, habitat %in% "offshore")) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  #scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
  #                  values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))
ggplot(data = subset(ssid.pca, habitat %in% "deep")) +
  theme_bw()+
  theme(panel.grid = element_blank()) +
  xlab(axis.labels[1]) +
  ylab(axis.labels[2]) +
  geom_point(aes(MDS1, MDS2, fill = admix), shape = 21, size = 5, alpha = 0.8)+
  #scale_fill_manual(name = "Ecomorph", labels = c("Nearshore", "Offshore", "Deep 1", "Deep 2"), 
  #                 values = c('palegreen1','plum1','skyblue','dodgerblue4')) +
  labs(fill = "Ecomorph", size = 20) +
  theme(axis.title = element_text(size =16),axis.text = element_text(size = 14),
        legend.title = element_text(size=16), legend.text = element_text(size = 14))


