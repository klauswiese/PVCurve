#Curva Presión - Volumen
#K. Wiese Julio 2018
###################################################################

PVCurveMean <- function(x,y){
  
  #Craer un data frame
  df <- data.frame(CRA = x, IPH = y, RSQ = rep(NA, length(y)))
  
  for (i in length(df$IPH):1){
    df$RSQ[i] <- cor(df$IPH[length(df$IPH):i], df$CRA[length(df$IPH):i])
  }
  
  Inflexion <- which.max(abs(df$RSQ[1:(length(df$RSQ)-2)]))
  Punto <- c(df$CRA[Inflexion], df$IPH[Inflexion])
  Punto2 <- c(df$CRA[Inflexion - 1], df$IPH[Inflexion - 1])
  
  PuntoMedio <- c((df$CRA[Inflexion] + df$CRA[Inflexion - 1])/2, 
                  (df$IPH[Inflexion] + df$IPH[Inflexion - 1])/2)
  
  
  
  tabla <- data.frame(rbind(df[Inflexion:dim(df)[1],1:2], PuntoMedio))
  
  
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
  points(PuntoMedio[1], PuntoMedio[2], pch=20, col="red")
  text(PuntoMedio[1] + 0.05, PuntoMedio[2] + 0.05, expression(paste(psi, rho, "= 0")))
  
  ####################################################
  #Inverso del potencial osmótico en plena turgencia #
  ####################################################
  IPOPT <- lm2$coefficients[1]
  points(0, IPOPT, pch=10)
  text(0.01, IPOPT + 0.1, expression(paste("1 / ", psi^100)))
  
  ####################################################
  #Inverso del potencial osmótico a turgencia cero   #
  ####################################################
  IPOTC <- PuntoMedio[2]
  points(0, PuntoMedio[2] , pch=8)
  text(0 +0.1, IPOTC, expression(paste("1/", psi^0)))
  
  #############################################
  #Contenido hídrico relativo a turgencia cero#
  #############################################
  CHRTC <- PuntoMedio[1] 
  points(CHRTC, 0 , pch=6)
  text(CHRTC, 0.15, expression(CHR^0))
  
  #############################
  #Volumen hídrico simplástico#
  #############################
  VHS <- (0 - lm2$coefficients[1]) / lm2$coefficients[2]
  points(VHS, 0 , pch=4)
  text(VHS, 0.15, expression(V_s))
  

  Resultados <- list(InversoPotencialOsmoticoPlenaTurgencia = as.numeric(IPOPT), 
                     InversoPotencialOsmoticoTurgenciaCero=as.numeric(IPOTC), 
                     ContenidoHidricoRelativoTurgenciaCero=as.numeric(CHRTC), 
                     VolumenHidricoSimplastico=as.numeric(VHS))
  return(Resultados)
  print(Resultados)
  
}

