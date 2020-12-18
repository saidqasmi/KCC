## Script permettant definir l'environnement spline avec n noeuds.
## Les notations sont celles de paper.ps, le 09/03/09

# Inputs : need a weighting function w to compute locally adaptative splines 


base_spline = function(n,w=1,rho=1){

## Base spline : 
x = 1:n	# Les knots
rn = n+2;	# Dimension de la base spline
X =matrix(data = 0,nrow = rn,ncol = 1) 
## Il faut chercher une facon de decrire la base qui ne pose pas (trop) 
##de probl??mes num??riques ; pour cela j'utilise des fonctions pr??-existantes de scilab.
##Pour d??crire une spline, j'utilise ses valeurs aux knots (n coefs), et ses valeurs au milieu des premier et dernier intervalles (X(2) et X($-1)), soit (n+2) coordonn??es.
X[1,1] = x[1]
X[2,1] = mean(x[1:2])
X[3:(rn-2),1] = x[2:(n-1)]
X[rn-1,1] = mean(x[(n-1):n])
X[rn,1] = x[n]

Y0 =diag(x=1,nrow = rn,ncol = rn)
Y1 = Y2 = Y3 = matrix(data=0,nrow=rn,ncol=rn)
for (i in 1:rn){
f <- splinefun(X,Y0[,i],method = "fmm")## Calcul de la fonction spline 	
## Calcul des valeurs (f0), dérivée (f1), dérivée seconde (f2) et troisième (f3) aux points Xi.
Y1[,i] = f(X,deriv=1);
Y2[,i] = f(X,deriv=2);
Y3[,i] = f(X,deriv=3);
}


## Checks on w
if (length(w)==1) { w=rep(1,n) }
nw = length(w)
if (nw != n){ stop("Length of w() is not correct") }
# Construct the final weights ww
w_interv = apply(cbind(w[1:(n-1)],w[2:n]),1,mean)
ww = w_interv[c(1,1:(n-1),n-1)]


## Matrice B :
  Zn = Y0[c(1, 3:(rn-2), rn),]
## Matrice G :
  Gn = G0=matrix(data = 0,nrow=dim(t(Zn)%*%Zn)[1],ncol = dim(t(Zn)%*%Zn)[2]) 
             
             interv = X[2:rn] - X[1:(rn-1)]
             interv =matrix(data = interv,nrow = length(interv),ncol = 1)
             Interv = interv%*% matrix(data = 1, nrow=1, ncol=length(Y2[1,]) )  
             
for(s1 in 1:rn){
      #nombre de noeuds copies en colonnes des deriv??e secondes et troisi??mes 
      #en les noeuds sauf le dernier pour la s1 ieme fct de base spline
      f2 = Y2[1:(dim(Y2)[1]-1),s1]%*% matrix(data = 1, nrow =1,ncol =length(Y2[1,])  )  
      f3 = Y3[1:(dim(Y3)[1]-1),s1]%*% matrix(data = 1, nrow =1,ncol =length(Y2[1,])  )  
      #deriv?? des fonction de base arrang??e en colonnes 
      #sur tout les noeud sauf le dernier
      g2 = Y2[1:(dim(Y2)[1]-1),]
      g3 = Y3[1:(dim(Y3)[1]-1),]
      #voire notes
      Int = f2 * g2 *Interv + (f2 * g3 + g2 * f3)*Interv^2/2 + f3 * g3 *Interv^3/3	## Pour comprendre ce calcul, il faut revenir ?? la d??finition de G et d??velopper le produit des d??riv??es secondes sur chaque intervalle...
     #sommation sur les colonnes de Int
       Gn[s1,] =apply(ww*Int, 2, sum)   
}

Zn <<- Zn
Gn <<- Gn
#Hn = t(Zn)%*%Zn + rho * Gn
#Gamma_n <<- Zn %*% solve(Hn) %*% t(Zn)
}



