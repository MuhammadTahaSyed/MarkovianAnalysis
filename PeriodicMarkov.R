library(markovchain)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(gridExtra)

setwd('/home/taha/Documents/MarkovChains/Steady') # input

SeasonalMarkov <- function(filename, acc = 3) {
  zz <- read.csv(filename)
  datez <- as.Date(zz$datex)
  montx <- month(datez)
  zza <- data.frame(zz, montx)
  
  states <- c("Dry", "Normal", "Wet")
  tarray <- array(0, dim = c(3, 3, 12))
  
  for (i in 1:12) {
    k <- ifelse(i == 12, 1, i + 1)
    from <- as.numeric(unlist(subset(zza, montx == i, select = spistatets)))
    to <- as.numeric(unlist(subset(zza, montx == k, select = spistatets)))
    if (length(from) > length(to)) from <- from[-1]
    if (length(to) > length(from)) to <- to[-1]
    mat <- table(from, to)
    normmat <- sweep(mat, 1, rowSums(mat), "/")
    tarray[,,i] <- normmat
  }
  
  mclist <- lapply(1:12, function(i) {
    mat <- prop.table(tarray[,,i], 1)
    new("markovchain", states = states, transitionMatrix = mat, name = paste("Month", month.abb[i]))
  })
  names(mclist) <- month.abb
  
  results <- lapply(1:12, function(i) {
    P <- diag(length(states))
    for (j in 0:11) {
      idx <- (i + j - 1) %% 12 + 1
      P <- P %*% mclist[[idx]]@transitionMatrix
    }
    mc <- new("markovchain", states = states, transitionMatrix = P, name = paste("Annual Chain =", month.name[i]))
    ss <- steadyStates(mc)
    mrt <- 1 / ss
    mfpt <- round(meanFirstPassageTime(mc,'Normal'), 3)
    list(Stationary = setNames(as.numeric(ss), states),
         MRT = setNames(as.numeric(mrt), states),
         MFPT = mfpt)
  })
  names(results) <- month.name
  return(results)
}

PSI = SeasonalMarkov('markstates_PRISM_SSMI.csv')
CSI = SeasonalMarkov('markstates_CNRM_SSMI.csv')
USI = SeasonalMarkov('markstates_UKESM_SSMI.csv')

levels = c(
  "October", "November", "December",
  "January", "February", "March", "April", "May", "June",
  "July", "August", "September"
)

# Convert Stationary values into a data frame
PSIdf <- do.call(rbind, lapply(names(PSI), function(m) {
  data.frame(Month = m, State = names(PSI[[m]]$Stationary), Value = PSI[[m]]$Stationary)
}))
PSIdf$Month <- factor(PSIdf$Month,levels)

CSIdf <- do.call(rbind, lapply(names(CSI), function(m) {
  data.frame(Month = m, State = names(CSI[[m]]$Stationary), Value = CSI[[m]]$Stationary)
}))
CSIdf$Month <- factor(CSIdf$Month,levels)

USIdf <- do.call(rbind, lapply(names(USI), function(m) {
  data.frame(Month = m, State = names(USI[[m]]$Stationary), Value = USI[[m]]$Stationary)
}))
USIdf$Month <- factor(USIdf$Month,levels)


#png("PRSIM  SSMI Seasonal Probabilities.png", width = 1800, height = 600, res = 100)

p1 <- ggplot(PSIdf, aes(x = Month, y = State, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = round(Value, 3)), size = 3) +
  scale_fill_gradient(low = "white", high = "cyan") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PRISM SSMI Probabilities", fill = "Probability")

p2 <- ggplot(transform(PSIdf, Value = as.numeric(as.character(Value)) -
                         as.numeric(as.character(CSIdf$Value))),
             aes(x = Month, y = State, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = round(Value, 3)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PRISM & CNRM Comparision", fill = "Difference")

p3 <- ggplot(transform(PSIdf, Value = as.numeric(as.character(Value)) -
                         as.numeric(as.character(USIdf$Value))),
             aes(x = Month, y = State, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = round(Value, 3)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "PRISM & UKESM Comparision", fill = "Difference")

grid.arrange(p2, p1, p3,ncol = 3)
#dev.off()


Pmrt <- do.call(rbind, lapply(names(PSI), function(m) {
  data.frame(Month = m, State = names(PSI[[m]]$MRT), Value = PSI[[m]]$MRT)
}))

Pmrt$Month <- factor(Pmrt$Month, levels)

Cmrt <- do.call(rbind, lapply(names(CSI), function(m) {
  data.frame(Month = m, State = names(CSI[[m]]$MRT), Value = CSI[[m]]$MRT)
}))

Cmrt$Month <- factor(Cmrt$Month, levels)

Umrt <- do.call(rbind, lapply(names(USI), function(m) {
  data.frame(Month = m, State = names(USI[[m]]$MRT), Value = USI[[m]]$MRT)
}))

Umrt$Month <- factor(Umrt$Month, levels)

#png("PRSIM  SSMI Seasonal MRT.png", width = 1800, height = 600, res = 100)

m1 = ggplot(Pmrt, aes(x = Month, y = State, fill = Value)) +
    geom_tile() +
    geom_text(aes(label = round(Value, 3)), size = 3) +
    scale_fill_gradient(low = "white", high = "lightgreen") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PRISM SSMI Mean Recurrence Time", fill = "MRT")

m2 = ggplot(transform(Pmrt, Value = as.numeric(as.character(Value)) -
                        as.numeric(as.character(Cmrt$Value))),
            aes(x = Month, y = State, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = round(Value, 3)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Difference PRISM vs CNRM", fill = "MRT")

m3= ggplot(transform(Pmrt, Value = as.numeric(as.character(Value)) -
                           as.numeric(as.character(Umrt$Value))),
               aes(x = Month, y = State, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = round(Value, 3)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Difference PRISM vs UKESM", fill = "MRT")

grid.arrange(m2,m1,m3, ncol =3)
#dev.off()




# Extract MFPT values into tidy data frame
Pmfpt <- do.call(rbind, lapply(names(PSI), function(m) {
  data.frame(
    Month = m,
    From = c("Dry", "Wet"),
    MFPT = c(PSI[[m]]$MFPT["Dry"], PSI[[m]]$MFPT["Wet"])
  )
}))
Cmfpt <- do.call(rbind, lapply(names(CSI), function(m) {
  data.frame(
    Month = m,
    From = c("Dry", "Wet"),
    MFPT = c(CSI[[m]]$MFPT["Dry"], CSI[[m]]$MFPT["Wet"])
  )
}))
Umfpt <- do.call(rbind, lapply(names(USI), function(m) {
  data.frame(
    Month = m,
    From = c("Dry", "Wet"),
    MFPT = c(USI[[m]]$MFPT["Dry"], USI[[m]]$MFPT["Wet"])
  )
}))

# Ensure month order
Pmfpt$Month <- factor(Pmfpt$Month, levels)

#png("PRSIM  SSMI MFPT.png", width = 1800, height = 600, res = 100)
# Plot heatmap
fp1 = ggplot(Pmfpt, aes(x = Month, y = From, fill = MFPT)) +
    geom_tile(color = "white") + 
    geom_text(aes(label = round(MFPT, 3)), color = "black", size = 3) +
    scale_fill_gradient(low = "white", high = "#20B2AA") +
    labs(title = "PRISM SSMI Mean First Passage Time to Normal State", x = "Month", y = "From State") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

fp2 = ggplot(transform(Pmfpt, MFPT = MFPT - Cmfpt$MFPT), aes(x = Month, y = From, fill = MFPT)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = round(MFPT, 3)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "PRISM vs Downscaled CNRM", x = "Month", y = "From State") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fp3 = ggplot(transform(Pmfpt, MFPT = MFPT - Umfpt$MFPT), aes(x = Month, y = From, fill = MFPT)) +
  geom_tile(color = "white") + 
  geom_text(aes(label = round(MFPT, 3)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "PRISM vs Downscaled UKESM", x = "Month", y = "From State", fill = 'Delta MFPT') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(fp2,fp1,fp3, ncol =3)
#dev.off()


