# Load required packages
library(markovchain)
library(ggplot2)
library(dplyr)

# Set working directory
path <- "/home/taha/Documents/MarkovChains/Seasonal" # input directory
setwd(path)
out <- "/home/taha/Documents/MarkovChains/Expanded"


# Function to expand and compute metrics for any file
Expand <- function(file) {
  tarray <- readRDS(file)
  states <- c("Dry", "Normal", "Wet")
  months <- month.abb
  n <- length(states) * length(months)
  expanded_states <- paste0(rep(states, each = 12), "_", rep(months, times = 3))
  mat <- matrix(0, n, n, dimnames = list(expanded_states, expanded_states))
  
  for (m in 1:12) {
    next_m <- ifelse(m == 12, 1, m + 1)
    P <- sweep(tarray[,,m], 1, rowSums(tarray[,,m]), "/")
    for (i in 1:3) for (j in 1:3) {
      r <- (i - 1) * 12 + m
      c <- (j - 1) * 12 + next_m
      mat[r, c] <- P[i, j]
    }
  }
  
  mc <- new("markovchain", states = expanded_states, transitionMatrix = mat)
  mfpt <- round(meanFirstPassageTime(mc), 3)
  ss <- round(as.numeric(steadyStates(mc)), 3)
  names(ss) <- states(mc)
  rec <- round(1 / ss, 3)
  
  list(mfpt = mfpt, steady = ss, recurrence = rec)
}

# Apply function to all models
PRISM <- Expand("PSSMI_tarray.rds")
CNRM <- Expand("CNRM_tarray.rds")
UKESM <- Expand("UKESM_tarray.rds")

# Combine recurrence data
df <- rbind(
  data.frame(State = names(PRISM$recurrence), Recurrence = PRISM$recurrence,
             Category = sub("_.*", "", names(PRISM$recurrence)), Model = "PRISM"),
  data.frame(State = names(CNRM$recurrence), Recurrence = CNRM$recurrence,
             Category = sub("_.*", "", names(CNRM$recurrence)), Model = "CNRM"),
  data.frame(State = names(UKESM$recurrence), Recurrence = UKESM$recurrence,
             Category = sub("_.*", "", names(UKESM$recurrence)), Model = "UKESM")
)

# Order states and months
month_order <- c("Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep") # Seasonal Sort
state_order <- c(paste0("Dry_", month.abb), paste0("Normal_", month.abb), paste0("Wet_", month.abb))
df$State <- factor(df$State, levels = state_order)

# Separate plots for Dry, Normal, Wet
categories <- c("Dry","Normal","Wet")

setwd(out) # output directory

for (cat in categories) {
  df_cat <- df %>%
    filter(Category == cat) %>%
    mutate(Month = sub(".*_", "", State),
           Month = factor(Month, levels = month_order))
  
  # Uncomment below to save PNG
  png(paste0("Recurrence_", cat, ".png"), width = 1800, height = 900, res = 150)
  print(
    ggplot(df_cat, aes(x = Month, y = Recurrence, fill = Model)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 1) +
      scale_fill_manual(values = c("PRISM" = "#35cc42", "CNRM" = "#6d77f7", "UKESM" = "#ff6045"),labels = c(
        "PRISM" = "PRISM", "CNRM" = "Downscaled CNRM","UKESM" = "Downscaled UKESM")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste("Recurrence Times for", cat, "States"),
           x = "Month", y = "Recurrence Time (months)")
  )
  dev.off()
}

# Steady-state probabilities
df2 <- rbind(
  data.frame(State = names(PRISM$steady), Recurrence = PRISM$steady,
             Category = sub("_.*", "", names(PRISM$steady)), Model = "PRISM"),
  data.frame(State = names(CNRM$steady), Recurrence = CNRM$steady,
             Category = sub("_.*", "", names(CNRM$steady)), Model = "CNRM"),
  data.frame(State = names(UKESM$steady), Recurrence = UKESM$steady,
             Category = sub("_.*", "", names(UKESM$steady)), Model = "UKESM")
)

for (cat in categories) {
  df_cat2 <- df2 %>%
    filter(Category == cat) %>%
    mutate(Month = sub(".*_", "", State),
           Month = factor(Month, levels = month_order))
  
  # Uncomment below to save PNG
  png(paste0("Steady_", cat, ".png"), width = 1800, height = 900, res = 150)
  print(
    ggplot(df_cat2, aes(x = Month, y = Recurrence, fill = Model)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 1) +
      scale_fill_manual(values = c("PRISM" = "#35cc42", "CNRM" = "#6d77f7", "UKESM" = "#ff6045"), labels = c(
        "PRISM" = "PRISM", "CNRM" = "Downscaled CNRM","UKESM" = "Downscaled UKESM")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste("Steady-State Probabilities for", cat, "States"),
           x = "Month", y = "Probability")
  )
  dev.off()
}

# MFPT extraction
extract_mfpt <- function(mfpt_matrix, model_name) {
  months <- month.abb
  data <- data.frame()
  for (m in months) {
    for (cat in c("Dry", "Wet")) {
      from_state <- paste0(cat, "_", m)
      to_state <- paste0("Normal_", m)
      mfpt_val <- mfpt_matrix[from_state, to_state]
      data <- rbind(data, data.frame(Month = m, Category = cat, MFPT = mfpt_val, Model = model_name))
    }
  }
  return(data)
}

df_mfpt <- rbind(
  extract_mfpt(PRISM$mfpt, "PRISM"),
  extract_mfpt(CNRM$mfpt, "CNRM"),
  extract_mfpt(UKESM$mfpt, "UKESM")
)
df_mfpt$Month <- factor(df_mfpt$Month, levels = month_order)

for (cat in c("Dry","Wet")) {
  df_cat_mfpt <- df_mfpt %>%
    filter(Category == cat)
  
  # Uncomment below to save PNG
  png(paste0("MFPT_", cat, ".png"), width = 1800, height = 900, res = 150)
  print(
    ggplot(df_cat_mfpt, aes(x = Month, y = MFPT, fill = Model)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
      geom_text(aes(label = round(MFPT, 3)), 
                position = position_dodge(width = 0.8), 
                angle = 90, vjust = 0.5, hjust = 2, size = 3, color = "white") +
      scale_fill_manual(values = c("PRISM" = "#35cc42", "CNRM" = "#6d77f7", "UKESM" = "#ff6045" ), labels = c(
        "PRISM" = "PRISM", "CNRM" = "Downscaled CNRM","UKESM" = "Downscaled UKESM")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste("Mean First Passage Time to Normal State from", cat, "State"),
           x = "Month", y = "Mean First Passage Time (months)")
  )
  dev.off()
}