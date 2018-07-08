#Curva Presión - Volumen
#K. Wiese Julio 2018
###################################################################

PVCurvelm <- function(x,y){
  
  #Crear un data frame
  df <- data.frame(CRA = x, IPH = y, RSQ = rep(NA, length(y)))
  
  #Calcular la relacion entre puntos usando R cuadrado de atras hacia adelante
  for (i in length(df$IPH):1){
    df$RSQ[i] <- cor(df$IPH[length(df$IPH):i], df$CRA[length(df$IPH):i])
  }
  
  #Encuentra el valor m'aximo de R cuadrado a partir del cual caen los valores de R cuadrado
  Inflexion <- which.max(abs(df$RSQ[1:(length(df$RSQ)-2)]))
  Punto <- c(df$CRA[Inflexion], df$IPH[Inflexion])#Punto maximo
  Punto2 <- c(df$CRA[Inflexion - 1], df$IPH[Inflexion - 1])#Punto de caida con respecto al maximo
  
  #Regresion entre los puntos en donde ocurre la caida
  t0 <- df[(Inflexion-4):length(df$CRA),]
  
  #t0 <-data.frame(rbind(Punto, Punto2))#data frame con los puntos donde ocurre la caida
  #names(t0) <- c("CRA", "IPH")#Nombres variables
  
  #lm0 <-lm(IPH ~ CRA, data=t0)#Regresion
  
  rhs <- function(x, b0, b1) {
    b0 + x^b1
  }
  
  lmpower <- nls(IPH ~ rhs(CRA, intercept, power), data = t0 ,
             start = list(intercept = 0, power = -2), trace = T)
  
  #Grafico de modelo 
  #plot(t0$CRA, predict(lmpower,t0$CRA),lty = 1, col = "blue", type="l")
  #lines(t0$CRA, predict(lmpower,t0$CRA),lty = 1, col = "blue", type="l")
  #Evaluacion de 10 puntos en la funcion generada a partir de los valores de la caida
  CRAlm <- seq(Punto2[1], Punto[1], length = 5)
  IPHlm <- predict(lmpower, list(CRA = CRAlm))
  Nuevos <- data.frame(cbind(CRAlm, IPHlm))
  row.names(Nuevos) <- paste("P", 1:5, sep="")
  names(Nuevos) <- c("CRA", "IPH")
  #Creacion de tabla con nuevos valores
  dfN <- rbind(Nuevos, df[Inflexion:length(df$CRA),1:2])
  
  
    #Calcular la relacion entre puntos usando R cuadrado de atras hacia adelante
  for (i in length(dfN$IPH):1){
    dfN$RSQ[i] <- cor(dfN$IPH[length(dfN$IPH):i], dfN$CRA[length(dfN$IPH):i])
  }
  
  
  #Encuentra el valor m'aximo de R cuadrado a partir del cual caen los valores de R cuadrado
  Inflexion2 <- which.max(abs(dfN$RSQ[1:(length(dfN$RSQ)-2)]))
  
  #plot(dfN$CRA, dfN$IPH)
  #lines(t0$CRA, predict(lmpower,t0$CRA), lty = 1, col = "blue")
  
  tabla <- dfN[Inflexion2:dim(dfN)[1],1:2]
  
  #Curva 2
  lm2 <- lm(IPH ~ CRA, data = tabla)
  
  #Distribución de Puntos
  plot(x, y, xlab="1 - CRA", ylab=expression(paste("1 / ", psi)), pch=16, 
       main="Curva Presión - Volumen \n Punto Medio",
       xlim=c(-0.001, max(x)+0.5), 
       ylim=c(-0.001, max(y)))
  abline(v=0)
  abline(h=0)
  grid(col="black")
  box()
  #abline(lm1, col="red")
  abline(lm2, col="darkgreen")
  #legend("topright", c("Curva 1", "Curva 2"), col=c("red", "darkgreen"), lty=1)
  #Índices derivados de la curva presión - volumen
  
  
  ####################################################
  #Inverso del potencial osmótico en plena turgencia #
  ####################################################
  PI <- c()
  points(tabla$CRA[1], tabla$IPH[1], pch=20, col="red")
  text(tabla$CRA[1] + 0.05, tabla$IPH[1] + 0.05, expression(paste(psi, rho, "= 0")))
  
  ####################################################
  #Inverso del potencial osmótico en plena turgencia #
  ####################################################
  IPOPT <- lm2$coefficients[1]
  points(0, IPOPT, pch=10)
  text(0.01, IPOPT + 0.1, expression(paste("1 / ", psi^100)))
  
  ####################################################
  #Inverso del potencial osmótico a turgencia cero   #
  ####################################################
  IPOTC <- predict(lm2)[1]
  points(0, IPOTC , pch=8)
  text(0 +0.1, IPOTC, expression(paste("1/", psi^0)))
  
  #############################################
  #Contenido hídrico relativo a turgencia cero#
  #############################################
  CHRTC <- tabla$CRA[1] 
  points(CHRTC, 0 , pch=6)
  text(CHRTC, 0.15, expression(CHR^0))
  
  #############################
  #Volumen hídrico simplástico#
  #############################
  VHS <- (0 - lm2$coefficients[1]) / lm2$coefficients[2]
  points(VHS, 0 , pch=4)
  text(VHS, 0.15, expression(V_s))
  
  ######################
  #Potencial de presión#   
  ######################
  #PP <- y[2] - (coef(lm2)[2] * x[2] + coef(lm2)[1])
  lmpower2 <- nls(IPH ~ rhs(CRA, intercept, power), data = df ,
                 start = list(intercept = 0, power = -2), trace = T)
  lines(df$CRA, predict(lmpower2,df$CRA),lty = 1, col = "blue", type="l")
  
  #Ecuación de la recta 2
  #Y <- coef(lm2)[2] * x + coef(lm2)[1]
  
  #legend("bottomright", c(expression(paste("1 / ", psi^100)), 
  #                        expression(paste("1/", psi^0)), 
  #                        expression(CHR^0), 
  #                        expression(V_s)), 
  #       pch=c(10,8,6,4))
  
  Resultados <- list(InversoPotencialOsmoticoPlenaTurgencia = as.numeric(IPOPT), 
                     InversoPotencialOsmoticoTurgenciaCero=as.numeric(IPOTC), 
                     ContenidoHidricoRelativoTurgenciaCero=as.numeric(CHRTC), 
                     VolumenHidricoSimplastico=as.numeric(VHS)) 
                     #PotencialPresion=as.numeric(PP))
  return(Resultados)
  print(Resultados)  
}

