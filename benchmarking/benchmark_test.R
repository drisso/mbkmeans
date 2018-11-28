# package
library(SingleCellExperiment)
library(mbkmeans)
library(DelayedMatrixStats)
library(ClusterR)
library(rhdf5)
library(TENxPBMCData)
library(scater)
library(ggplot2)


# Speed compare
mat <- HDF5Array("benchmarking/pbmc3k_rectangular.h5", name = "counts")
mat_t<-t(mat)
data_hdf5 <- mat[1:2000,1:1000]
data_matrix <- as.matrix(data_hdf5)
data_sce <- SingleCellExperiment(assays = list(counts = data_hdf5))

#kmeans++

system.time(clusterR_res <- MiniBatchKmeans(t(data_inmem), clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))

system.time(mbk_res<-mbkmeans(data_matrix,clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))

system.time(mbk_res<-mbkmeans(data_hdf5,clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))

system.time(mbk_res<-mbkmeans(data_sce, reduceMethod = NA, clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10))

#random
system.time(clusterR_res <- MiniBatchKmeans(t(data_inmem), clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10,initializer = "random"))

system.time(mbk_res<-mbkmeans(data_inmem,clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10,initializer = "random"))

system.time(mbk_res<-mbkmeans(data_hdf5,clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10,initializer = "random"))

system.time(mbk_res<-mbkmeans(data_sce, reduceMethod = NA, clusters = 3, batch_size = 100, init_fraction = .1, max_iters = 10,initializer = "random"))




#Transposition Burden
mat <- HDF5Array("benchmarking/pbmc3k_rectangular.h5", name = "counts")
mat<-as.matrix(mat)
mat_t<-t(mat)
sample_data<-as.matrix(mat[1:2000,1:1000])

#memory for the full matrix

#HTml

hdf_matrix<-as(sample_data,"HDF5Matrix")
hdf_matrix_t<-as(t(sample_data),"HDF5Matrix")

#mini_bathch function

mini_bat<-numeric(100)
for(i in 1:100){
    mini_bat[i]<-system.time(mini_batch(hdf_matrix_t,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mini_bat_mean<-mean(mini_bat)
mini_bat_sd<-sd(mini_bat)


mini_bat_t<-numeric(100)
for(i in 1:100){
    mini_bat_t[i]<-system.time(mini_batch(t(hdf_matrix),clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mini_bat_t_mean<-mean(mini_bat_t)
mini_bat_t_sd<-sd(mini_bat_t)

#mbkmeans function

mbk<-numeric(100)
for(i in 1:100){
    mbk[i]<-system.time(mbkmeans(hdf_matrix_t,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mbk_mean<-mean(mbk)
mbk_sd<-sd(mbk)


mbk_t<-numeric(100)
for(i in 1:100){
    mbk_t[i]<-system.time(mbkmeans(t(hdf_matrix_t),clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mbk_t_mean<-mean(mbk_t)
mbk_t_mean_sd<-sd(mbk_t)


boxplot(mini_bat,mini_bat_t,mbk,mbk_t,
        main ="Boxplot of Transposition Burden Time",
        xlab = "Type",
        ylab = "Time")
axis(side = 1,at=seq(1,4),label = c("mini_bat","mini_bat_t","mbk","mbk_t"))


trans_mean<-rbind(mini_bat_mean,mini_bat_t_mean,mbk_mean,mbk_t_mean)
trans_sd<-rbind(mini_bat_sd,mini_bat_t_sd,mbk_sd,mbk_t_mean_sd)
trans_total<-data.frame(cbind(trans_mean,trans_sd))
colnames(trans_total)<-c("Mean","SD")

ggplot(trans_total, aes(x=rownames(trans_total), y=Mean)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD)) +
    labs(title = "Mean and SD of total time of functions with/without transposition",
          x = "Funciton",
          y = "Time")



#whole matrix

mini_bat_whole<-numeric(100)
for(i in 1:100){
    mini_bat_whole[i]<-system.time(mini_batch(mat,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mini_bat_whole_mean<-mean(mini_bat_whole)
mini_bat_whole_sd<-sd(mini_bat_whole)


mini_bat_whole_t<-numeric(100)
for(i in 1:100){
    mini_bat_whole_t[i]<-system.time(mini_batch(t(mat),clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mini_bat_whole_t_mean<-mean(mini_bat_whole_t)
mini_bat_whole_t_sd<-sd(mini_bat_whole_t)

#mbkmeans function

mbk_whole<-numeric(100)
for(i in 1:100){
    mbk_whole[i]<-system.time(mbkmeans(mat_t,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mbk_whole_mean<-mean(mbk_whole)
mbk_whole_sd<-sd(mbk_whole)



mbk_whole_t<-numeric(100)
for(i in 1:100){
    mbk_whole_t[i]<-system.time(mbkmeans(t(mat),clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))[3]
    i = i+1
}
mbk_whole_t_mean<-mean(mbk_whole_t)
mbk_whole_t_sd<-sd(mbk_whole_t)



boxplot(mini_bat_whole,mini_bat_whole_t,mbk_whole,mbk_whole_t,
        main ="Boxplot of Transposition Burden Time(Whole data)",
        xlab = "Type",
        ylab = "Time")
axis(side = 1,at=seq(1,4),label = c("mini_bat_whole","mini_bat_whole_t","mbk_whole","mbk_whole_t"))


trans_mean_whole<-rbind(mini_bat_whole_mean,mini_bat_whole_t_mean,
                        mbk_whole_mean,mbk_whole_t_mean)
trans_sd_whole<-rbind(mini_bat_whole_sd,mini_bat_whole_t_sd,
                      mbk_whole_sd,mbk_whole_t_sd)
trans_total_whole<-data.frame(cbind(trans_mean_whole,trans_sd_whole))
colnames(trans_total_whole)<-c("Mean","SD")

ggplot(trans_total_whole, aes(x=rownames(trans_total_whole), y=Mean)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD)) +
    labs(title = "Mean and SD of total time of functions with/without transposition(Whole data)",
         x = "Funciton",
         y = "Time")







#memory  --change packages and funciotn
# gc(mini_bat<-mini_batch(hdf_matrix_t,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))
#
# gc(mini_bat<-mini_batch(t(hdf_matrix),clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))

#mbkmeans function

# #memory
# gc(mini_bat<-mbkmeans(hdf_matrix,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))
#
# gc(mini_bat<-mbkmeans(t(hdf_matrix_t),clusters = 3,batch_size = 100, init_fraction = .1, max_iters = 10))



#For memory, the transposition of data won't lead to any change. However, the used time is different for different funciton. For mini_batch funtion, the time for two methods is almost same. However, for the mbkmeans function, the calculation of simple data spends less time compared with the calculation of transposition.


#Default values


#max_iters (ClusterR=100,  mbkmeans = 10)
#
# time_max<-numeric(10)
# for(i in 1:10){
# time_max[i]<-system.time(mbkmeans(hdf_matrix,clusters = 3,batch_size = 100, init_fraction = .1, max_iters = i))[3]
# }
#
# plot(1:10,time_max,type = "b",xlab = "max_iters",ylab = "Time")


#init_fraction (ClusterR -1, mbkmeans =0.25, suggestion = 0.8)
#
# time_init<-numeric(10)
# for(i in 1:10){
# time_init[i]<-system.time(mbkmeans(hdf_matrix,clusters = 3,batch_size = 100, init_fraction =i/10, max_iters = 10))[3]
# }
#
# plot(seq(0.1,1,0.1),time_init,type = "b",xlab = "init_fraction",ylab = "Time")



