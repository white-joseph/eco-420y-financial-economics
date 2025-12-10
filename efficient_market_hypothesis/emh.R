############################################################
## ECO 420Y Project – Testing EMH for Bitcoin (Part II)
## R code for Sections 1–4
############################################################

## Packages
library(quantmod)
library(tseries)
library(xtable)

## Paths
base_dir <- "C:/Repositories/eco-420y-financial-economics/efficient_market_hypothesis"
plot_dir <- file.path(base_dir, "plots")

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## Section 1: Download Bitcoin price data from FRED
############################################################

sym <- "CBBTCUSD"
getSymbols(sym, src = "FRED")

btc_price_xts <- get(sym)
colnames(btc_price_xts) <- "Price"
btc_price_xts <- na.omit(btc_price_xts)

head(btc_price_xts)

############################################################
## Section 2: Statistical analysis
## - daily returns
## - plots
## - descriptive statistics
## - unit root tests
############################################################

## Daily log returns
btc_ret_xts <- diff(log(btc_price_xts$Price))
btc_ret_xts <- na.omit(btc_ret_xts)
colnames(btc_ret_xts) <- "Return"

btc_price <- as.numeric(btc_price_xts$Price)
btc_ret   <- as.numeric(btc_ret_xts$Return)

## Plot price and return (saved to file)
png(file = file.path(plot_dir, "bitcoin_price_return.png"),
    width = 1200, height = 800)
par(mfrow = c(2, 1))

plot(btc_price_xts,
     main = "Bitcoin Daily Price (USD)",
     ylab = "Price",
     xlab = "Date")

plot(btc_ret_xts,
     main = "Bitcoin Daily Log Return",
     ylab = "Return",
     xlab = "Date")

par(mfrow = c(1, 1))
dev.off()

## Descriptive statistics
descriptive_stats <- function(x) {
  x <- x[is.finite(x)]
  c(
    mean   = mean(x),
    sd     = sd(x),
    min    = min(x),
    q25    = quantile(x, 0.25),
    median = median(x),
    q75    = quantile(x, 0.75),
    max    = max(x)
  )
}

price_stats  <- descriptive_stats(btc_price)
return_stats <- descriptive_stats(btc_ret)

stats_table <- rbind(Price = price_stats,
                     Return = return_stats)
round(stats_table, 4)

## LaTeX table for descriptive statistics
latex_tab <- xtable(
  stats_table,
  caption = "Descriptive Statistics for Bitcoin Price and Daily Log Return",
  label   = "tab:btc_stats",
  digits  = 4
)

tex_path <- file.path(base_dir, "btc_descriptive_stats.tex")

# Capture the LaTeX table as a character vector
latex_lines <- capture.output(
  print(latex_tab, type = "latex", include.rownames = TRUE)
)

# Write those lines to the .tex file (same style as your TWFE tables)
writeLines(latex_lines, tex_path, useBytes = TRUE)


## ADF tests for price and return
adf_price  <- adf.test(btc_price, alternative = "stationary", k = 1)
adf_price

adf_return <- adf.test(btc_ret, alternative = "stationary", k = 1)
adf_return

############################################################
## Section 3: Rolling AR(1) t-values (window size = 28)
## y_t = phi_0 + phi_1 y_{t-1} + error_t
############################################################

r     <- btc_ret
n     <- length(r)
rlag1 <- c(NA, r[1:(n - 1)])

d <- data.frame(r = r, rlag1 = rlag1)
d <- na.omit(d)
d$tr <- 1:nrow(d)

window <- 28
Td     <- nrow(d)

t_values <- numeric(Td - window + 1)

j <- 1
while (j + window - 1 <= Td) {
  subd <- subset(d, tr >= j & tr <= j + window - 1)
  m    <- lm(r ~ rlag1, data = subd)
  t_values[j] <- summary(m)$coef["rlag1", "t value"]
  j <- j + 1
}

## Save rolling AR(1) t-values plot
png(file = file.path(plot_dir, "bitcoin_rolling_tvalues.png"),
    width = 1000, height = 700)

plot(t_values, type = "o",
     xlab = "Window index",
     ylab = "t value of phi_1",
     main = "Rolling AR(1) t-values (Bitcoin returns, window = 28)")
abline(h = -1.96, lty = 2)
abline(h =  1.96, lty = 2)

dev.off()

############################################################
## Section 4: Rolling Box–Ljung p-values (window size = 28)
############################################################

r_full   <- btc_ret
n_full   <- length(r_full)
n_win    <- n_full - window + 1
p_values <- numeric(n_win)

j <- 1
while (j + window - 1 <= n_full) {
  subr <- r_full[j:(j + window - 1)]
  bj   <- Box.test(subr, lag = 1, type = "Ljung-Box")
  p_values[j] <- bj$p.value
  j <- 1 + j
}

## Save rolling Box–Ljung p-values plot
png(file = file.path(plot_dir, "bitcoin_rolling_pvalues.png"),
    width = 1000, height = 700)

plot(p_values, type = "o",
     xlab = "Window index",
     ylab = "p-value (Box–Ljung, lag 1)",
     main = "Rolling Box–Ljung p-values (Bitcoin returns, window = 28)",
     ylim = c(0, 1))
abline(h = 0.05, lty = 2)

dev.off()