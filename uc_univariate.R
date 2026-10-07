## =====================================================================
##  uc_univariate.R
##  Trend-cycle decomposition of US real GDP (FRED: GDPC1)
##     y_t       = mu_t + c_t
##     mu_{t+1}  = delta_t + mu_t + eta_t,              eta_t ~ N(0, s2_eta)   [I(1) trend]
##     c_{t+1}   = phi1 c_t + phi2 c_{t-1} + kappa_t,  kappa_t ~ N(0, s2_kappa) [AR(2) cycle]
##  State: alpha_t = (mu_t, c_t, c_{t-1})'
##  (the "UC-0" model of Morley, Nelson & Zivot 2003: uncorrelated trend and cycle shocks)
## =====================================================================
source("ssm_functions.R"); source("plot_helpers.R")
data_dir <- "../data"; fig_dir <- "../figures"
covid_missing <- c("2020Q2", "2020Q3")     # treated as missing observations

## ---- data
gdp <- read_fred_quarterly(file.path(data_dir, "GDPC1.csv"))
gdp <- gdp[gdp$quarter >= "1948Q1", ]
y   <- 100 * log(gdp$value)
y[gdp$quarter %in% covid_missing] <- NA
tt  <- quarter_to_time(gdp$quarter); n <- length(y)
cat("Sample:", gdp$quarter[1], "-", gdp$quarter[n], " n =", n, "\n")

## ---- model
## Drift with a break in 1973Q1 (productivity slowdown; Perron & Wada 2009):
##   delta_t = delta_0 + delta_1 * 1{t >= 1973Q1}.  The trend stays I(1).
break_dates <- "1973Q1"
make_D <- function(breaks) { D <- matrix(1, n, 1); for (b in breaks) D <- cbind(D, as.numeric(gdp$quarter >= b)); D }

build_uc <- function(th, D, y1) {
  th <- unname(th); k <- ncol(D)
  phi <- ar2_from_pacf(0.995 * tanh(th[k + 3]), 0.995 * tanh(th[k + 4]))
  Tm  <- matrix(c(1, 0,      0,
                  0, phi[1], phi[2],
                  0, 1,      0), 3, 3, byrow = TRUE)
  R   <- matrix(c(1, 0,
                  0, 1,
                  0, 0), 3, 2, byrow = TRUE)
  Q   <- diag(exp(th[k + 1:2]))
  P1  <- matrix(0, 3, 3); P1[1, 1] <- 1e7                                   # diffuse trend
  P1[2:3, 2:3] <- stationary_P(Tm[2:3, 2:3], (R %*% Q %*% t(R))[2:3, 2:3])  # stationary cycle
  list(Z = matrix(c(1, 1, 0), 1, 3), H = matrix(0), T = Tm, R = R, Q = Q,
       c = cbind(D %*% th[1:k], 0, 0),                                      # time-varying drift
       a1 = c(y1, 0, 0), P1 = P1,
       par = c(setNames(th[1:k], paste0("delta", 0:(k - 1))),
               s2_eta = exp(th[k + 1]), s2_kappa = exp(th[k + 2]), phi1 = phi[1], phi2 = phi[2]))
}

fit_uc <- function(D) {
  k <- ncol(D)
  negll <- function(th) {
    v <- try(-loglik(y, build_uc(th, D, y[1]), drop = 1), silent = TRUE)
    if (inherits(v, "try-error") || !is.finite(v)) 1e10 else v
  }
  best <- NULL                                   # multistart BFGS
  for (r1 in c(0.8, 0.95)) for (s2 in c(0.1, 0.5)) {
    th0 <- c(0.8, rep(-0.2, k - 1), log(s2), log(0.5), atanh(r1), atanh(-0.4))
    o <- optim(th0, negll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10))
    if (is.null(best) || o$value < best$value) best <- o
  }
  mod <- build_uc(best$par, D, y[1])
  ## standard errors of the natural parameters: inverse Hessian + delta method
  V_th <- pinv(num_hessian(negll, best$par))
  g <- function(th) build_uc(th, D, y[1])$par
  J <- sapply(seq_along(best$par), function(j) {
    e <- replace(rep(0, length(best$par)), j, 1e-6); (g(best$par + e) - g(best$par - e)) / 2e-6 })
  list(mod = mod, logL = -best$value, par = mod$par, se = sqrt(pmax(diag(J %*% V_th %*% t(J)), 0)))
}

fit0 <- fit_uc(make_D(character(0)))            # constant drift (for comparison)
fit  <- fit_uc(make_D(break_dates))              # drift break in 1973Q1
mod <- fit$mod; par_hat <- fit$par; se <- fit$se; logL <- fit$logL
phi <- par_hat[c("phi1", "phi2")]
roots <- polyroot(c(1, -phi[1], -phi[2]))
period <- if (phi[1]^2 + 4 * phi[2] < 0) 2 * pi / acos(phi[1] / (2 * sqrt(-phi[2]))) else NA
print(round(cbind(estimate = par_hat, se = se), 4))
cat("logL =", logL, " (constant drift:", fit0$logL, ", LR =", 2 * (logL - fit0$logL), ")\n")
cat("inverse AR roots (moduli):", 1 / Mod(roots), "  period (q):", period, "\n")

## ---- filter / smoother
kf <- kalman_filter(y, mod); sm <- state_smoother(y, mod, kf)
cyc_s  <- sm$ahat[, 2]; cyc_sd <- sqrt(pmax(sm$V[2, 2, ], 0))
cyc_f  <- kf$att[, 2];  cyc_fsd <- sqrt(pmax(kf$Ptt[2, 2, ], 0))
trend_s <- sm$ahat[, 1]

kf0 <- kalman_filter(y, fit0$mod); cyc0 <- state_smoother(y, fit0$mod, kf0)$ahat[, 2]

## ---- diagnostics on standardised innovations (t = 2..n)
e  <- std_innov(kf)[-1]; e <- e[!is.na(e)]
lb8  <- Box.test(e, lag = 8,  type = "Ljung-Box", fitdf = length(par_hat) - 1)
lb12 <- Box.test(e, lag = 12, type = "Ljung-Box", fitdf = length(par_hat) - 1)
sk <- mean((e - mean(e))^3) / sd(e)^3; ku <- mean((e - mean(e))^4) / var(e)^2
JB <- length(e) / 6 * (sk^2 + (ku - 3)^2 / 4)
h  <- floor(length(e) / 3); Hh <- sum(tail(e, h)^2) / sum(head(e, h)^2)
diag_tab <- c(Q8 = lb8$statistic, p8 = lb8$p.value, Q12 = lb12$statistic, p12 = lb12$p.value,
              JB = JB, pJB = 1 - pchisq(JB, 2), skew = sk, kurt = ku, H = Hh, h = h,
              pH = 2 * min(pf(Hh, h, h), 1 - pf(Hh, h, h)))
print(round(diag_tab, 3))

## ---- forecasts
hmax <- 16
fc <- ssm_forecast(mod, kf, hmax)

saveRDS(list(break_dates = break_dates, logL0 = fit0$logL, par0 = fit0$par, cyc0 = cyc0, quarter = gdp$quarter, tt = tt, y = y, par = par_hat, se = se, logL = logL, period = period,
             trend = trend_s, trend_sd = sqrt(pmax(sm$V[1, 1, ], 0)), cyc_s = cyc_s, cyc_sd = cyc_sd,
             cyc_f = cyc_f, cyc_fsd = cyc_fsd, diag = diag_tab, fc = fc, e = e),
        "results_univariate.rds")

## ---- figures ------------------------------------------------------
yr <- range(tt)
open_pdf(file.path(fig_dir, "fig_gdp_trend.pdf"), 10, 3.2)
par(mfrow = c(1, 2), cex = 0.9)
plot(tt, y, type = "n", xlab = "", ylab = "", main = "100 x log real GDP and smoothed trend")
shade_recessions(par("usr")[3:4])
lines(tt, y, col = UMgrey, lwd = 2.2); lines(tt, trend_s, col = UMdblue, lwd = 1.2)
legend("topleft", c(expression(y[t]), expression(hat(mu)[t])), col = c(UMgrey, UMdblue),
       lwd = c(2.2, 1.2), bty = "n", cex = 0.85)
dy <- c(NA, diff(y)) * 4
plot(tt, dy, type = "h", col = UMgrey, xlab = "", ylab = "", ylim = c(-12, 12),
     main = "Quarterly growth (annualised %), COVID quarters missing")
lines(tt, 4 * drop(make_D(break_dates) %*% par_hat[grep("delta", names(par_hat))]), col = UMorange, lwd = 1.8, lty = 2)
legend("bottomleft", c(expression(4 * Delta * y[t]), expression(4 * hat(delta)[t])),
       col = c(UMgrey, UMorange), lty = c(1, 2), lwd = c(1, 1.5), bty = "n", cex = 0.85)
invisible(dev.off())

open_pdf(file.path(fig_dir, "fig_gdp_cycle.pdf"), 10, 3.4)
yl <- range(c(cyc_s - 1.645 * cyc_sd, cyc_s + 1.645 * cyc_sd, cyc_f), na.rm = TRUE)
plot(tt, cyc_s, type = "n", ylim = yl, xlab = "", ylab = "% of trend",
     main = "Output gap (cycle): smoothed with 90% band, and filtered")
shade_recessions(yl)
band(tt, cyc_s - 1.645 * cyc_sd, cyc_s + 1.645 * cyc_sd, UMdblue)
abline(h = 0, col = UMgrey, lwd = 0.7)
lines(tt, cyc_f, col = UMorange, lwd = 1.1, lty = 2)
lines(tt, cyc_s, col = UMdblue, lwd = 1.6)
legend("bottomleft", c(expression(hat(c)[t] == E(c[t] * "|" * Y[n])), expression(a[t * "|" * t] == E(c[t] * "|" * Y[t]))),
       col = c(UMdblue, UMorange), lty = c(1, 2), lwd = c(1.6, 1.1), bty = "n", cex = 0.85, horiz = TRUE)
invisible(dev.off())

open_pdf(file.path(fig_dir, "fig_gdp_diag.pdf"), 11, 2.8)
par(mfrow = c(1, 3), cex = 0.85)
te <- tt[-1][!is.na(std_innov(kf)[-1])]
plot(te, e, type = "l", col = UMdblue, xlab = "", ylab = "", main = expression("standardised innovations " * v[t] / sqrt(F[t])))
abline(h = 0, col = UMgrey)
ac <- acf(e, lag.max = 16, plot = FALSE)$acf[-1]
plot(1:16, ac, type = "h", lwd = 6, lend = 1, col = UMdblue, ylim = c(-0.25, 0.25), xlab = "lag", ylab = "", main = "ACF")
abline(h = c(-1, 1) * 1.96 / sqrt(length(e)), col = UMorange, lty = 2); abline(h = 0)
hist(e, breaks = 30, freq = FALSE, col = alpha_col(UMdblue, 0.5), border = "white", main = "histogram vs N(0,1)", xlab = "", ylab = "")
curve(dnorm(x), add = TRUE, col = UMorange, lwd = 1.5)
invisible(dev.off())

open_pdf(file.path(fig_dir, "fig_gdp_forecast.pdf"), 10, 3.1)
sel <- tt >= 2012
tf <- tail(tt, 1) + (1:hmax) / 4
yl <- range(c(y[sel], fc$mean - 1.645 * fc$sd, fc$mean + 1.645 * fc$sd), na.rm = TRUE)
plot(tt[sel], y[sel], type = "n", xlim = c(2012, max(tf)), ylim = yl, xlab = "", ylab = "",
     main = "Forecasts: the Kalman filter without updates (90% intervals)")
band(tf, fc$mean - 1.645 * fc$sd, fc$mean + 1.645 * fc$sd, UMorange, 0.2)
lines(tt[sel], y[sel], col = UMgrey, lwd = 2.2); lines(tt[sel], trend_s[sel], col = UMdblue, lwd = 1.1)
lines(c(tail(tt, 1), tf), c(tail(y, 1), fc$mean), col = UMorange, lwd = 1.6)
legend("topleft", c(expression(y[t]), expression(hat(mu)[t]), expression(hat(y)[n + h * "|" * n])),
       col = c(UMgrey, UMdblue, UMorange), lwd = c(2.2, 1.1, 1.6), bty = "n", cex = 0.85)
invisible(dev.off())

## ---- HP filter as a Kalman smoother (integrated random walk + noise, q = 1/1600)
hp_mod <- list(Z = matrix(c(1, 0), 1, 2), H = matrix(1), T = matrix(c(1, 0, 1, 1), 2, 2),
               R = matrix(c(0, 1), 2, 1), Q = matrix(1 / 1600), a1 = c(y[1], 0), P1 = diag(c(1e7, 1e7)))
hp_kf <- kalman_filter(y, hp_mod); hp_trend <- state_smoother(y, hp_mod, hp_kf)$ahat[, 1]
## check against the penalised least-squares HP formula on the complete-data version
ycomp <- y; miss <- is.na(ycomp)
D <- diff(diag(n), differences = 2)
W <- diag(as.numeric(!miss)); ycomp[miss] <- 0
hp_pls <- solve(W + 1600 * crossprod(D), W %*% ycomp)
cat("max |KF smoother - HP penalised LS| =", max(abs(hp_trend - hp_pls)), "\n")
hp_cycle <- y - hp_trend
open_pdf(file.path(fig_dir, "fig_gdp_hp.pdf"), 10, 3.0)
yl <- range(c(cyc_s, hp_cycle), na.rm = TRUE)
plot(tt, cyc_s, type = "n", ylim = yl, xlab = "", ylab = "",
     main = sprintf("UC cycle vs HP cycle (correlation %.2f)", cor(cyc_s, hp_cycle, use = "complete.obs")))
shade_recessions(yl); abline(h = 0, col = UMgrey, lwd = 0.7)
lines(tt, hp_cycle, col = UMorange, lwd = 1.1, lty = 2); lines(tt, cyc_s, col = UMdblue, lwd = 1.5)
legend("bottomleft", c("UC: I(1) trend + AR(2) cycle", expression("HP, " * lambda == 1600)),
       col = c(UMdblue, UMorange), lty = c(1, 2), lwd = c(1.5, 1.1), bty = "n", cex = 0.85)
invisible(dev.off())

## ---- Kalman filter at work on GDP: estimates and conditional variances of the cycle
covid <- which(gdp$quarter %in% covid_missing)
open_pdf(file.path(fig_dir, "fig_gdp_kf.pdf"), 9, 2.8)
par(mfrow = c(1, 2), cex = 0.9)
sel <- tt >= 2005
yl <- range(c(kf$a[1:n, 2][sel], cyc_f[sel], cyc_s[sel]))
plot(tt[sel], cyc_s[sel], type = "n", ylim = yl, xlab = "", ylab = "% of trend", main = "Cycle: predicted, filtered, smoothed")
shade_recessions(yl); abline(h = 0, col = UMgrey, lwd = 0.7)
lines(tt[sel], kf$a[1:n, 2][sel], col = UMgrey, lwd = 1.2, lty = 3)
lines(tt[sel], cyc_f[sel], col = UMorange, lwd = 1.3)
lines(tt[sel], cyc_s[sel], col = UMdblue, lwd = 1.6)
legend("bottomleft", c(expression("predicted " * a[t]), expression("filtered " * a[t * "|" * t]), expression("smoothed " * hat(c)[t])),
       col = c(UMgrey, UMorange, UMdblue), lty = c(3, 1, 1), lwd = c(1.2, 1.3, 1.6), bty = "n", cex = 0.8)
Pp <- kf$P[2, 2, 1:n]; Pf <- kf$Ptt[2, 2, ]; Vs <- sm$V[2, 2, ]
plot(tt, Pp, type = "l", col = UMgrey, lwd = 1.3, lty = 3, ylim = c(0, max(Pp[-(1:4)]) * 1.05), xlab = "", ylab = "",
     main = "Conditional variances of the cycle")
lines(tt, Pf, col = UMorange, lwd = 1.3); lines(tt, Vs, col = UMdblue, lwd = 1.5)
legend("topleft", c(expression(P[t] == Var(c[t] * "|" * Y[t - 1])), expression(P[t * "|" * t] == Var(c[t] * "|" * Y[t])),
       expression(V[t] == Var(c[t] * "|" * Y[n]))), col = c(UMgrey, UMorange, UMdblue), lty = c(3, 1, 1), lwd = 1.4, bty = "n", cex = 0.8)
arrows(2016, 5.6, 2020, max(Pp[covid]) + 0.1, length = 0.05, col = UMorange); text(2016, 5.6, "COVID quarters missing", pos = 2, cex = 0.75, col = UMorange)
invisible(dev.off())
cat("cycle variances: P steady ~", round(median(Pp), 3), " Ptt ~", round(median(Pf), 3), " V ~", round(median(Vs), 3),
    " | at COVID P:", round(Pp[covid], 2), " V:", round(Vs[covid], 2), "\n")

## ---- missing observations: the COVID quarters
open_pdf(file.path(fig_dir, "fig_gdp_missing.pdf"), 9, 2.7)
par(mfrow = c(1, 2), cex = 0.9)
sel <- tt >= 2016 & tt < 2024
yhat <- kf$a[1:n, 1] + kf$a[1:n, 2]                       # one-step prediction Z a_t
trend_sd <- sqrt(pmax(sm$V[1, 1, ], 0))
gdp_all <- 100 * log(gdp$value)
yl <- range(c(gdp_all[sel], trend_s[sel] + 1.645 * trend_sd[sel], trend_s[sel] - 1.645 * trend_sd[sel]))
plot(tt[sel], y[sel], type = "n", ylim = yl, xlab = "", ylab = "", main = expression("Level: " * y[t] * ", prediction " * Z * a[t] * ", smoothed trend"))
rect(tt[covid[1]] - 0.1, yl[1] - 5, tt[covid[2]] + 0.1, yl[2] + 5, col = alpha_col(UMorange, 0.12), border = NA)
band(tt[sel], (trend_s - 1.645 * trend_sd)[sel], (trend_s + 1.645 * trend_sd)[sel], UMdblue, 0.18)
points(tt[covid], gdp_all[covid], pch = 4, col = UMorange)
lines(tt[sel], trend_s[sel], col = UMdblue, lwd = 1.4)
lines(tt[sel], yhat[sel], col = UMgrey, lty = 2)
points(tt[sel], y[sel], pch = 16, cex = 0.55, col = "black")
legend("topleft", c(expression(y[t]), "excluded (COVID)", expression(Z * a[t]), expression(hat(mu)[t] * " (90%)")),
       pch = c(16, 4, NA, NA), lty = c(NA, NA, 2, 1), col = c("black", UMorange, UMgrey, UMdblue), bty = "n", cex = 0.75)
yl <- range(c(cyc_s[sel] + 1.645 * cyc_sd[sel], cyc_s[sel] - 1.645 * cyc_sd[sel]))
plot(tt[sel], cyc_s[sel], type = "n", ylim = yl, xlab = "", ylab = "% of trend", main = "Smoothed cycle with 90% band")
rect(tt[covid[1]] - 0.1, yl[1] - 5, tt[covid[2]] + 0.1, yl[2] + 5, col = alpha_col(UMorange, 0.12), border = NA)
band(tt[sel], (cyc_s - 1.645 * cyc_sd)[sel], (cyc_s + 1.645 * cyc_sd)[sel], UMdblue, 0.18)
abline(h = 0, col = UMgrey, lwd = 0.7); lines(tt[sel], cyc_s[sel], col = UMdblue, lwd = 1.5)
invisible(dev.off())

open_pdf(file.path(fig_dir, "fig_gdp_nobreak.pdf"), 10, 3.0)
yl <- range(c(cyc0, cyc_s), na.rm = TRUE)
plot(tt, cyc0, type = "n", ylim = yl, xlab = "", ylab = "% of trend", main = "Smoothed cycle: constant drift vs drift break in 1973Q1")
shade_recessions(yl); abline(h = 0, col = UMgrey, lwd = 0.7); abline(v = 1973, col = UMorange, lty = 3)
lines(tt, cyc0, col = UMorange, lwd = 1.3, lty = 2); lines(tt, cyc_s, col = UMdblue, lwd = 1.5)
legend("bottomleft", c("constant drift", "drift break 1973Q1"), col = c(UMorange, UMdblue), lty = c(2, 1),
       lwd = c(1.3, 1.5), bty = "n", cex = 0.85, horiz = TRUE)
invisible(dev.off())
cat("done\n")
