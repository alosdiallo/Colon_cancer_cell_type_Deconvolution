## generate figures for project

level5means <- colMeans(DHColon_HTD)
level4means <- colMeans(DHColon_L4)
level3means <- colMeans(DHColon_L3)
level2means <- colMeans(DHColon_L2)
level1means <- colMeans(DHColon_L1)

mean1 <- colMeans(df.cell %>% 
                    filter(MLH1==0,Distant_Mets==0) %>% 
                    select(2:16))
mean2 <- colMeans(df.cell %>% 
                    filter(MLH1==0,Distant_Mets==1) %>% 
                    select(2:16))
mean3 <- colMeans(df.cell %>% 
                    filter(MLH1==1,Distant_Mets==0) %>% 
                    select(2:16))
mean4 <- colMeans(df.cell %>% 
                    filter(MLH1==1,Distant_Mets==1) %>% 
                    select(2:16))
dhmc_combined <- rbind(mean1, mean2, mean3, mean4)
# _______________________________________________

m1 <- colMeans(tcga.data %>% 
                 filter(MLH1==0,Distant_Mets==0) %>% 
                 select(2:18))
m2 <- colMeans(tcga.data %>% 
                 filter(MLH1==0,Distant_Mets==1) %>% 
                 select(2:18))
m3 <- colMeans(tcga.data %>% 
                 filter(MLH1==1,Distant_Mets==0) %>% 
                 select(2:18))
m4 <- colMeans(tcga.data %>% 
                 filter(MLH1==1,Distant_Mets==1) %>% 
                 select(2:18))
dhmc_combined <- rbind(m1, m2, m3, m4)


# _____________________________________________
# EWAS Figs

DHMC_hypo <- res %>% 
  filter(logFC < -1) %>% 
  filter(P.Value < 0.05/50)

DHMC_hyper <- res %>% 
  filter(logFC > 1) %>% 
  filter(P.Value < 0.05/50)

DHMC_hypo_adj <- res2 %>% 
  filter(logFC < -1) %>% 
  filter(P.Value < 0.05/50)

DHMC_hyper_adj <- res %>% 
  filter(logFC > 1) %>% 
  filter(P.Value < 0.05/50)



