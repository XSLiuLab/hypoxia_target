library(dplyr)
###将 metalab 的格式转化为 xml 格式
###调用 matlab 进行并行
# devtools::install_github("muschellij2/matlabr")
library(matlabr)
options(matlab.path = "/home/data/sdb/wt/matlab2023/bin/")
# have_matlab()

all_files <- list.files("/home/data/sdb/wt/pan-cancer-GEMs/all_models/")
done <- list.files("/home/data/sdb/wt/GEMs_TCGA/")
remains <- all_files[which(!(all_files %in% paste0(gsub(".xml","",done),"_TINIT.mat")))]
remains <- paste0("/home/data/sdb/wt/pan-cancer-GEMs/all_models/",remains)
library(doParallel)
library(foreach)
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
res <- foreach(
  i = remains,
  .packages = "matlabr"
) %dopar% {
  code = c("addpath(genpath('/home/data/sdb/wt/RAVEN'));", 
           paste0("tmp = load('",i,"');"), 
           "exportModel(tmp.TINITmodel, sprintf('/home/data/sdb/wt/GEMs_TCGA/%s.xml', tmp.TINITmodel.id));"
  )
  tmp = run_matlab_code(code)
}
parallel::stopCluster(cl = my.cluster)

###将 xml 转化成酶网络
library(dplyr)
library(Met2Graph)
xml_files <- list.files("/home/data/sdb/wt/GEMs_TCGA/",full.names = F)
done <- list.files("/home/data/sdb/wt/MetGraph/EnzGraphs/",
                   pattern = "_enzymes_based_graph.ncol")
remains <- xml_files[which(!(xml_files %in% paste0(gsub("_enzymes_based_graph.ncol","",done),".xml")))]
remains <- paste0("/home/data/sdb/wt/GEMs_TCGA/",remains)
library(doParallel)
library(foreach)
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
res <- foreach(
  i = remains,
  .packages = "Met2Graph"
) %dopar% {
  Met2Graph::Met2EnzGraph(i, rmMets=TRUE, 
                          outDir="/home/data/sdb/wt/MetGraph/", 
                          outFormat="ncol")
}
parallel::stopCluster(cl = my.cluster)



