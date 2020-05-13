library(ggplot2)
library(openxlsx)
library(reshape2)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(Cairo)
library(igraph)

data <- read.table("E:\\thesis\\final_model_R2.txt",header=T)
all_phe <- read.xlsx("E:\\thesis\\H2.xlsx",)
all_phe <- all_phe[3:203,c(1:3,13:14)]
colnames(all_phe) <- all_phe[1,]
all_phe <- all_phe[-1,]
all_phe[,5] <- as.numeric(all_phe[,5])

for (i in 2:200) {
  if (is.na(all_phe[i,1])){
    all_phe[i,1] <- all_phe[i-1,1]
  }
}

exemplar_trait <- all_phe[which(all_phe$Heritability>0.35),]   ### 将遗传力大于0.35的性状作为研究对象
exemplar_trait <- exemplar_trait[,-c(2,5)]
colnames(exemplar_trait) <- c("Assay","phenotype","category")

#### data reshape
new_data <- melt(data,id.vars = c("phenotype","h2","SNP_num"),measure.vars = colnames(data)[4:27]) %>%
  merge(.,exemplar_trait,by="phenotype") 
 

temp<- data.frame()
for (i in 1:dim(new_data)[1]){
  if(str_detect(new_data[i,4],"lm")){
    a <- str_sub(new_data[i,4],1,-4)
    temp[i,1] <- "MLR"
  }else if(str_detect(new_data[i,4],"cv")){
    a <- str_sub(new_data[i,4],1,-4)
    temp[i,1] <- "ENR"
  }else if(str_detect(new_data[i,4],"svr")){
    a <- str_replace(new_data[i,4],"svr_","")
    temp[i,1] <- "SVR"
  }
  temp[i,2] <- a
}  ##### temp有两列，第一列为模型类型，第二列为数据类型

new_data <- new_data[,-4] 
new_data <- cbind(new_data[,1:4],temp,new_data[,5:6]) %>%
  na.omit(.)

for (i in 1:592){
  if(new_data[i,3]>=1500){
    new_data[i,9] <- "H"
  }else if(new_data[i,3]<=500){
    new_data[i,9] <- "L"
  }else if(new_data[i,3]>500 & new_data[i,3]<1500){
    new_data[i,9] <- "M"
  }
}
new_data <- new_data[,-3]

for (i in 1:592){
  if(new_data[i,2]>=0.6){
    new_data[i,9] <- "H"
  }else if(new_data[i,2]<=0.45){
    new_data[i,9] <- "L"
  }else if(new_data[i,2]>0.45 & new_data[i,2]<0.6){
    new_data[i,9] <- "M"
  }
}

colnames(new_data) <- c("phenotype","h2","prediction","model","data","assay","category","SNP_num","h2_class")

####################  针对不同的目的清洗出特定的数据框 并绘图 ####################

######    3.2	不同模型类型的比较    ######

test_num <- which(str_detect(new_data[,5],"test"))
photo_data1 <- new_data[-test_num,]  #### delete test data

photo_data1$data <- str_replace_all(photo_data1$data,"_train","")

p <- ggplot(photo_data1,aes(h2,prediction)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme(plot.margin = unit(c(5,0.3,5,0.25),"lines"),legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  facet_grid(~model,scales="fixed") +
  geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
  geom_point(size=1.5,aes(color=category)) +
  scale_colour_manual(values= c(behaviour="maroon",physiological="goldenrod",tissue="lightskyblue"))





#### 
test_num <- which(str_detect(new_data[,5],"test"))
photo_data1 <- new_data[-test_num,]  #### delete test data

photo_data1$data <- str_replace_all(photo_data1$data,"_train","")

 p <- ggplot(photo_data1,aes(h2,prediction)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme(plot.margin = unit(c(5,0.3,5,0.25),"lines"),legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  facet_grid(~model,scales="fixed") +
  geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
  geom_point(size=1.5,aes(color=category)) +
  scale_colour_manual(values= c(behaviour="maroon",physiological="goldenrod",tissue="lightskyblue"))

 p
 CairoPNG("E:\\thesis\\picture\\final\\paper\\3.22.png",width = 6.42,height = 4.68,units='in',dpi=600)
 plot(p)
 dev.off()
 
 

#### 特征变量数增加，两个线性模型的比较
data_num2 <- which(str_detect(new_data[,4],"SVR"))
photo_data2 <- new_data[-data_num2,]

photo_data2 <- within(photo_data2,SNP_num <- factor(SNP_num,levels = c("L","M","H")))  ###
with(photo_data2,levels(SNP_num))

p <- ggplot(photo_data2,aes(h2,prediction)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme(plot.margin = unit(c(5,0,5,0.25),"lines"),legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  facet_grid(model~SNP_num) +
  geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
  geom_point(size=1.5,aes(color=category)) +
  scale_colour_manual(values= c(behaviour="maroon",physiological="goldenrod",tissue="lightskyblue"))

CairoPNG("E:\\thesis\\picture\\final\\paper\\3.21.png",width = 6.61,height = 5.34,units='in',dpi=600)
plot(p)
dev.off()


########  不同数据类型之间的比较  #######

#### SNP模型
test_num <- which(str_detect(new_data[,5],"test"))
photo_data3 <- new_data[-test_num,]  #### delete test data

photo_data3$data <- str_replace_all(photo_data3$data,"_train","")

data_num3 <- which(str_detect(photo_data3[,5],"trait"))
photo_data3 <- photo_data3[-data_num3,]

photo_data3 <- within(photo_data3,data <- factor(data,levels = c("highSNP","randomSNP","random2SNP")))  ###
with(photo_data3,levels(data))

p <- ggplot(photo_data3,aes(h2,prediction)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme(plot.margin = unit(c(5,0.25,5,0.25),"lines"),legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  facet_grid(~data) +
  geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
  geom_point(size=1.5,aes(color=category)) +
  scale_colour_manual(values= c(behaviour="maroon",physiological="goldenrod",tissue="lightskyblue"))

p
CairoPNG("E:\\thesis\\picture\\final\\paper\\3.31.png",width = 6.34,height = 4.63,units='in',dpi=600)
plot(p)
dev.off()


####  trait模型
test_num <- which(str_detect(new_data[,5],"test"))  ###去掉test数据集
photo_data4 <- new_data[-test_num,]  #### delete test data

photo_data4$data <- str_replace_all(photo_data4$data,"_train","")

data_num4 <- which(photo_data4[,5]=="trait")
data_num5 <- which(photo_data4[,5]=="randomSNP")
photo_data4 <- photo_data4[c(data_num4,data_num5),]

photo_data4 <- within(photo_data4,data <- factor(data,levels = c("randomSNP","trait")))  ###
with(photo_data4,levels(data))

p1 <- ggplot(photo_data4,aes(h2,prediction)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  #theme(plot.margin = unit(c(0.8,1,8,0.5),"lines"),legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  theme(legend.position = "bottom") +
  facet_grid(~data) +
  geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
  geom_point(size=1.5,aes(color=category)) +
  scale_colour_manual(values= c(behaviour="maroon",physiological="goldenrod",tissue="lightskyblue"))

p2 <- ggplot(photo_data4,aes(data,prediction),color=data) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  scale_y_continuous(limits = c(0,0.75)) +
  xlab(NULL) +
  facet_grid(~model) +
  geom_boxplot(aes(fill=factor(data))) +
  scale_colour_manual(values = c(randomSNP="maroon",trait="goldenrod")) +
  labs(fill="data")

  
p <- multiplot(p2,p1,cols = 1)

CairoPNG("E:\\thesis\\picture\\final\\3.32.png",width = 4.95,height = 4.43,units='in',dpi=1200)
plot(p)
dev.off()



#### SNP + trait
test_num <- which(str_detect(new_data[,5],"test"))  ###去掉test数据集
photo_data5 <- new_data[-test_num,]  #### delete test data

photo_data5$data <- str_replace_all(photo_data5$data,"_train","")

data_num4 <- which(photo_data5[,5]=="highSNP_trait")
data_num5 <- which(photo_data5[,5]=="randomSNP_trait")
data_num6 <- which(photo_data5[,5]=="trait")
photo_data5 <- photo_data5[c(data_num4,data_num5,data_num6),]

photo_data5 <- within(photo_data5,data <- factor(data,levels = c("trait","highSNP_trait","randomSNP_trait")))  ###
with(photo_data5,levels(data))

p1 <- ggplot(photo_data5,aes(h2,prediction)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme(plot.margin = unit(c(1,0.25,1,0.25),"lines"),legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  facet_grid(~data) +
  geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
  geom_point(size=1.5,aes(color=category)) +
  scale_colour_manual(values= c(behaviour="maroon",physiological="goldenrod",tissue="lightskyblue"))

p2 <- ggplot(photo_data5,aes(data,prediction),color=data) +
  theme(legend.position = "bottom",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),axis.line = element_line(colour="black")) +
  scale_y_continuous(limits = c(0,1)) +
  xlab(NULL) +
  facet_grid(~model) +
  geom_boxplot(aes(fill=factor(data))) +
  scale_colour_manual(values = c(trait="maroon",highSNP_trait="goldenrod",randomSNP_trait="lightskyblue")) +
  labs(fill="data")

p <- multiplot(p2,p1,cols = 1)

CairoPNG("E:\\thesis\\picture\\final\\3.33.png",width = 5.94,height = 5.31,units='in',dpi=1200)
plot(p)
dev.off()




############  遗传力大小对预测的影响   #####
test_num <- which(str_detect(new_data[,5],"test"))
photo_data6 <- new_data[-test_num,]  #### delete test data
photo_data6$data <- str_replace_all(photo_data6$data,"_train","")


data_num5 <- which(photo_data6[,5]=="randomSNP_trait")
data_num6 <- which(photo_data6[,5]=="trait")
data_num7 <- which(photo_data6[,5]=="randomSNP")
photo_data6 <- photo_data6[c(data_num5,data_num6,data_num7),]

photo_data6 <- photo_data6[which(photo_data6$model=="SVR"),]

photo_data6 <- within(photo_data6,data <- factor(data,levels = c("randomSNP","trait","randomSNP_trait")))  ###
photo_data6 <- within(photo_data6,h2_class <- factor(SNP_num,levels = c("L","M","H")))
with(photo_data6,levels(data))

if(FALSE){
  p1 <- ggplot(photo_data6,aes(h2,prediction)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    theme(plot.margin = unit(c(1,0.25,1,0.25),"lines"),
          legend.position = "bottom") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "gray96",fill = "gray96"),
          axis.line = element_line(colour="black")) +
    theme(strip.text.y = element_text(size = 7) ) +
    facet_grid(~h2_class) +
    geom_abline(slope = 1,intercept = 0,colour="#990000",linetype="dashed") +
    geom_point(size=1.5,aes(color=data)) +
    scale_colour_manual(values= c(trait="maroon",highSNP_trait="goldenrod",randomSNP_trait="lightskyblue"))
  
} #不需要散点图


p2 <- ggplot(photo_data6,aes(data,prediction),color=data) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),
        axis.line = element_line(colour="black")) +
  scale_y_continuous(limits = c(0.25,0.85)) +
  xlab(NULL) +
  facet_grid(~h2_class) +
  geom_boxplot(aes(fill=factor(data))) +
  scale_colour_manual(values = c(randomSNP="yellow3",trait="maroon",randomSNP_trait="lightskyblue")) +
  labs(fill="data") +
  theme(legend.position = "bottom")


CairoPNG("E:\\thesis\\picture\\final\\3.41.png",width = 3.84,height = 2.31,units='in',dpi=600)
plot(p2)
dev.off()





###################  模型效果提升-堆叠柱状图   #############

data <- read.table("E:\\thesis\\final_model_R2.txt",header=T)
prediction_R2 <- data$randomSNP_trait_svr_train
h2 <- data$h2

p_data <- data.frame(value = c(prediction_R2,h2),
                     trait = c(rep(1:27,2)),
                     contion = c(rep("prediction_R2",27),rep("h2",27)))

p <- ggplot(p_data,aes(x=trait,y=value))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray96",fill = "gray96"),
        axis.line = element_line(colour="black")) +
  geom_bar(stat = "identity",position="dodge",aes(fill=contion),width = 0.5) +
  xlab("Trait") +
  ylab(NULL) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_x_continuous(breaks = seq(1,27,1))

p
CairoPNG("E:\\thesis\\picture\\final\\3.34.png",width = 7.33,height = 3.11,units='in',dpi=600)
plot(p)
dev.off()

################### 抽象模型网络图  ###############

trait_name <- str_c("T",1:10)
SNP_name   <- str_c("S",1:10)
E_name     <- str_c("E",1:10)
O_name     <- str_c("O",1:10)

edges_choose <- function(m){
  if (m=="T1"){
    class_num <- sample(3:4,4,replace = T)
  }else{
    class_num <- sample(1:2,4,replace = T)
  }
  cor_trait <- sample(trait_name,class_num[1])
  cor_SNP   <- sample(SNP_name,class_num[2])
  cor_E     <- sample(E_name,class_num[3])
  cor_O     <- sample(O_name,class_num[4])
  from      <- c(cor_trait,cor_SNP,cor_E,cor_O)
  temp      <- data.frame(from = from,
                          to   = rep(m,length(from)),
                          fre  = runif(length(from),min=0,max=0.4))
  return(temp)
}

edges <- data.frame()
for(i in c("T1",SNP_name,O_name)){
  data  <- edges_choose(i)
  edges <- rbind(edges,data)
}

from <- as.character(edges[,1])
to   <- as.character(edges[,2])
self_num <- which(from == to)
T1_num <- which(from=="T1")
edges <- edges[-unique(c(self_num,T1_num)),]  ### 去掉自连接的线以及从T1出发的线

for(i in 1:dim(edges)[1]){
  if(str_detect(edges[i,1],"T") & str_detect(edges[i,2],"T")){
    edges[i,4] <- 1
  }else if(TRUE %in% str_detect(edges[i,1],c("S","O","E")) & str_detect(edges[i,2],"T")){
    edges[i,4] <- 2
  }else if(str_detect(edges[i,1],c("T")) & TRUE %in% str_detect(edges[i,2],c("S","O"))){
    edges[i,4] <- 3
  }else{
    edges[i,4] <- 4
  }
  print(i)
} ### 对线进行分类，T to T1 is 1,SOE to T1 is 2,T to SO is 3, other are 4

colnames(edges)[4] <- "class"

vertices <- data.frame(
  name  = c(trait_name,SNP_name,E_name,O_name),
  class = c(1,rep(2,9),rep(3,10),rep(4,10),rep(5,10)), ##### 1,2,3,4分别对应性状、SNP、环境、其他
  size  = c(17,rep(14,9),runif(30,min=4,max=8))
)

graph_data <- graph_from_data_frame(edges,directed = T,vertices = vertices)

#set.seed(1)

V(graph_data)$size <- V(graph_data)$size
color1 <- c("firebrick","hotpink","lightskyblue","seagreen","mediumorchid")
V(graph_data)$color <- color1[V(graph_data)$class]

E(graph_data)$width <- E(graph_data)$fre
E(graph_data)$arrow.size=0.1
color2 <- c("royalblue1","springgreen2","darkorchid","gray")
E(graph_data)$color <- color2[E(graph_data)$class]
E(graph_data)$width <- E(graph_data)$fre *6

l <- layout.fruchterman.reingold(graph_data)
plot(graph_data,layout=l)

legend("bottomleft",
       legend = c("Aim","Trait","SNP","Environment","Others"),
       col = c("firebrick","hotpink","lightskyblue","seagreen","mediumorchid"),
       fill = c("firebrick","hotpink","lightskyblue","seagreen","mediumorchid"),
      #pch = 1,
       bty = "n")

CairoPNG("E:\\thesis\\picture\\final\\3.42.png",width = 10,height = 10,units='in',dpi=800)
plot(graph_data,layout=l)
dev.off()
