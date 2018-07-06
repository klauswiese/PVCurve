#Curva Presión - Volumen
#K. Wiese Julio 2018
###################################################################

PVCurve <- function(x,y){
  
  #Encontrar punto de inflexión
  Inflexion <- ede(x, y, 0)
  
  #Generar curvas según punto de inflexión
  #Curva 1
  lm1 <- lm(y[1:Inflexion[1]] ~ x[1:Inflexion[1]])
  
  #Curva 2
  lm2 <- lm(y[Inflexion[1]:length(y)] ~ x[Inflexion[1]:length(x)])
  
  
  #Distribución de Puntos
  plot(x, y, xlab="1 - CRA", ylab=expression(paste("1 / ", psi)), pch=16, main="Curva Presión - Volumen",
       xlim=c(-0.001, max(x)+0.5), ylim=c(-0.001, max(y)))
  abline(v=0)
  abline(h=0)
  grid(col="black")
  box()
  abline(lm1, col="red")
  abline(lm2, col="darkgreen")
  legend("topright", c("Curva 1", "Curva 2"), col=c("red", "darkgreen"), lty=1)
  #Índices derivados de la curva presión - volumen
  
  
  ####################################################
  #Inverso del potencial osmótico en plena turgencia #
  ####################################################
  PI <- c()
  points(x[Inflexion[1]], y[Inflexion[1]], pch=20, col="red")
  text(x[Inflexion[1]] + 0.05, y[Inflexion[1]] + 0.05, expression(paste(psi, rho, "= 0")))
  
  ####################################################
  #Inverso del potencial osmótico en plena turgencia #
  ####################################################
  IPOPT <- lm2$coefficients[1]
  points(0, IPOPT, pch=10)
  text(0.01, IPOPT + 0.1, expression(paste("1 / ", psi^100)))
  
  ####################################################
  #Inverso del potencial osmótico a turgencia cero   #
  ####################################################
  IPOTC <- predict(lm2, 1)[1]
  points(0, IPOTC , pch=8)
  text(0 +0.1, IPOTC, expression(paste("1/", psi^0)))
  
  #############################################
  #Contenido hídrico relativo a turgencia cero#
  #############################################
  CHRTC <- x[Inflexion[1]] 
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
  PP <- y[2] - (coef(lm2)[2] * x[2] + coef(lm2)[1])
  
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
                     VolumenHidricoSimplastico=as.numeric(VHS), 
                     PotencialPresion=as.numeric(PP))
  return(Resultados)
  print(Resultados)
  
}

