#Datos
curva <- read.csv("data/Results.csv")

#Librerias
#remotes::install_github("klauswiese/PVCurve")
library(remotes)
library(PVCurve)

#Inverso IPH
CRA <- 1 - curva$CRA
IPH <- 1 / curva$Psi

#Curva
PVCurve(CRA, IPH)

#Ajuste de la curva
m <- nls(IPH~I(CRA^power),start = list(power = -1), trace = T)

#Explorar sumario resultado
summary(m)


#Gráfico de curva ajustada vr curva conocida
power <- round(summary(m)$coefficients[1], 3)
power.se <- round(summary(m)$coefficients[2], 3)
plot(IPH ~ CRA, main = "Fitted power model", sub = "Blue: fit; green: known")
s <- seq(0, 1, length = 100)
lines(s, s^3, lty = 2, col = "green")
lines(s, predict(m, list(CRA = s)), lty = 1, col = "blue")
text(0, 0.5, paste("y =x^ (", power, " +/- ", power.se,
                          ")", sep = ""), pos = 4)


rhs <- function(x, b0, b1) {
     b0 + x^b1
}

m.2 <- nls(IPH ~ rhs(CRA, intercept, power), 
             start = list(intercept = 0, power = 0), trace = T)


plot(IPH ~ CRA, main = "Fitted power model, with intercept",
          sub = "Blue: fit; magenta: fit w/o intercept; green: known")
abline(h = 0, lty = 1, lwd = 0.5)
lines(s, s^3, lty = 2, col = "green")
lines(s, predict(m.2, list(CRA = s)), lty = 1, col = "blue")
lines(s, predict(m, list(CRA = s)), lty = 2, col = "magenta")
segments(CRA, IPH, CRA, fitted(m.2), lty = 2, col = "red")


AIC(m)
AIC(m.2)



#LOG
m.l <- lm(log(IPH) ~ CRA)
summary(m.l)


plot(log(IPH) ~ CRA, xlab = "x", ylab = "log(y+.1)", main = "Log-linear fit")
abline(m.l)
text(0, 0.4, pos = 4, paste("log(y) = ", round(coefficients(m.l)[1],3), "+", round(coefficients(m.l)[2], 3)))


par(mfrow = c(2, 2))
plot(m.l)
par(mfrow = c(1, 1))


#Exponential model
m.e <- nls(IPH ~ I(exp(1)^(a + b * CRA)), start = list(a = 0, b = 1), trace = T)
summary(m.e)


a <- round(summary(m.e)$coefficients[1, 1], 4)
b <- round(summary(m.e)$coefficients[2, 1], 4)
plot(IPH ~ CRA, main = "Fitted exponential function", sub = "Blue: fit; green: known")
s <- seq(0, 1, length = 100)
lines(s, s^3, lty = 2, col = "green")
lines(s, predict(m.e, list(CRA = s)), lty = 1, col = "blue")
text(0, 0.5, paste("y =e^ (", a, " + ", b, " * x)", sep = ""), pos = 4)
lines(s, predict(m.2, list(CRA = s)), lty = 1, col = "red")


#Craer un data frame
df <- data.frame(CRA = CRA, IPH = IPH, RSQ = rep(NA, length(IPH)))

for (i in length(df$IPH):1){
   df$RSQ[i] <- cor(df$IPH[length(df$IPH):i], df$CRA[length(df$IPH):i])
   plot(abs(df$RSQ))
   abline(h=max(abs(df$RSQ[1:(length(df$RSQ)-2)])), col="red")
}


