d = read.table("MD_cyclopropane.dat", header=T)
par(mgp=c(2,1,0))

anim = function() {
  for (i in seq(1,nrow(d),6)) {
    TITLE = sprintf("MD du cyclopropane\nt  =%8.4f ps,   Epot  = %8.4f kcal/mol",d$t[i], d$E_pot[i])
    plot(0, type="n", xlim=c(0.0,1.4), ylim=c(0.0,1.3), pch=19, cex=3.5,main=TITLE,xlab="x (Ångström)",ylab="y (Ångström)")
    
    # Nouvelles positions:
    points(c(0,d$xB[i],d$xC[i]),c(0,0,d$yC[i]), pch=19, cex=3.5)
    
    # Lines reliant les carbones:
    lines(c(0,d$xC[i]), c(0,d$yC[i]) ,lwd=2)
    lines(c(0,d$xB[i]), c(0,0) ,lwd=2)
    lines(c(d$xC[i], d$xB[i]), c(d$yC[i],0) ,lwd=2)
    
    # Position initiale des carbones:
    points(d$xC[1], d$yC[1], pch=3, col="grey", cex=2, lwd=2)
    points(d$xC[1], d$yC[1], pch=21, col="grey", cex=3, lwd=2)
    points(d$xB[1], 0, pch=3, col="grey", cex=2, lwd=2)
    points(d$xB[1], 0, pch=21, col="grey", cex=3, lwd=2)
    points(0, 0, pch=3, col="grey", cex=2, lwd=2)
    points(0, 0, pch=21,col="grey", cex=3, lwd=2)
    
    # Barycentre:
    points(d$xG[i], d$yG[i], pch=4, col="black", cex=2, lwd=2)
    
    points(d$xC[1]+0.08, d$yC[1]+0.08, pch="C", col="black", cex=1.5, lwd=2)
    points(d$xB[1]+0.08, 0.08, pch="B", col="black", cex=1.5, lwd=2)
    points(0.08, 0.08, pch="A", col="black", cex=1.5, lwd=2)
    
    #Traçage du déplacement:
    for (j in seq(1,i)) {
      points(c(d$xB[j],d$xC[j]),c(0,d$yC[j]), pch=".", cex=2, col="red")
    }
    
    # Forces:
    if (d$F_xB[i] != 0) arrows(d$xB[i],0,d$xB[i]+(d$F_xB[i]/150),0,col=2,lwd=2,length=0.1,angle=20)
    if (d$F_xC[i] != 0 && d$F_yC[i] != 0) arrows(d$xC[i],d$yC[i],d$xC[i]+(d$F_xC[i]/150),d$yC[i]+(d$F_yC[i]/150),col=2,lwd=2,length=0.1,angle=20)
  }
}
saveGIF(expr = anim(), interval=.1, outdir = getwd(), movie.name="anim_MD.gif")

png("Energies_MD.png")
inter = 0:500
plot(d$t[inter],d$E_tot[inter],type="l",xlab="t (ps)",ylab="Energie (kcal/mol)",col="black", ylim=c(0,max(d$E_tot)))
points(d$t[inter],d$E_pot[inter],col="red",type="l")
points(d$t[inter],d$E_cin[inter],col="green",type="l")
legend("topleft",max(E_tot),c("Etot","Epot","Ekin"),col=1:3,lty=1)
dev.off()
