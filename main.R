rm(list=ls())
set.seed(11)



library(microbenchmark)
library(ggplot2)
library(gtools)



#fichiers pour génerer des HMM 
source(file = "continuous_hmm_generator.R")
source(file = "discrete_hmm_generator.R")


#fichiers pour la premiere problematique
source(file = "forward_discrete.R")
source(file = "forward_continuous.R")

source(file = "backward_discrete.R")
source(file = "backward_continuous.R")

source(file = "naif_discrete.R")


#fichiers pour la seconde problematique 
source(file = "em_discrete.R")
source(file = "em_continuous.R")


#fichier pour la troisieme problematique
source(file = "gibbs_sampler.R")


#fichier pour la quatrieme problematique
source(file = "viterbi_discrete.R")
source(file = "naif_viterbi_discrete.R")




#Parametres modele - Exemple 1
mu <- rep(.5, 2)
a <- matrix(c(1/3, 2/3, .5, .5), nrow = 2, ncol = 2, byrow = T)
p <- matrix(c(1/4, 3/4, .5, .5), nrow = 2, ncol = 2, byrow = T)
states <- 0:1


#Parametres modele - EM vs Gibbs

mu_gibbs <- c(.2, .6, .2)
p_gibbs <- matrix(c(.6, .3, .1, .1, .8, .1, .1, .3, .6), nrow = 3, ncol = 3, byrow = T)
m_gibbs <- c(-2, 0, 2)
sigma_gibbs <- 1


#Generation HMM
hmm_discrete_8 <- discrete_hmm_generator(mu, p, a, 8)
hmm_discrete_1000 <- discrete_hmm_generator(mu, p, a, 1000)

hmm_continuous_8 <- continuous_hmm_generator(mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs, 8)
hmm_continuous_1000 <- continuous_hmm_generator(mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs, 1000)

#-----------------------------------------------------

#Problematique 1: Test algorithme forward et backward 

#cadre discret

forward_discrete(hmm_discrete_8$observed, mu, p, a)$proba
forward_discrete(hmm_discrete_1000$observed, mu, p, a)$proba

backward_discrete(hmm_discrete_8$observed, mu, p, a)$proba
backward_discrete(hmm_discrete_1000$observed, mu, p, a)$proba


#cadre continu
forward_continuous(hmm_continuous_8$observed, mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs)$proba
forward_continuous(hmm_continuous_1000$observed, mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs)$proba

backward_continuous(hmm_continuous_8$observed, mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs)$proba
backward_continuous(hmm_continuous_1000$observed, mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs)$proba


#On retrouve bien les meme resultats pour la probabilite d apparition d une
#chaine d observation en utilisant les algorithmes forward et backward.

#Comparaison complexite des algorithme forward, backward et naif (Figure 1)


#premiere courbe
c_for <- rep(0, 8)
c_back <- rep(0, 8)
c_naif <- rep(0, 8)

for(i in 2:9){
  c_for[i - 1] <- median(microbenchmark(forward_discrete(hmm_discrete_8$observed[1 : i], mu, p, a)[[1]])$time)
  c_back[i - 1] <- median(microbenchmark(backward_discrete(hmm_discrete_8$observed[1 : i], mu, p, a)[[1]])$time)
  c_naif[i - 1] <- median(microbenchmark(naif_discrete(hmm_discrete_8$observed[1 : i], mu, p, a, states))$time)
}
abscisse<- 2 : 9

df <- data.frame(x = abscisse, y1 = c_for*10**-6, y2 = c_back*10**-6, y3 = c_naif*10**-6)

ggplot(df, aes(x)) + geom_line(aes(y = y1, colour = "forward"), linewidth = 2) + geom_line(aes(y = y2, colour = "backward"), linewidth = 2, linetype = "dashed") + geom_line(aes(y = y3, colour = "naïf"), linewidth = 2) + scale_x_continuous(name = "Longueur de la chaîne") + scale_y_continuous(name = "Temps d'éxectution (ms)") + theme_bw() +
  scale_color_manual(name="", values = c("forward" = "brown1", "backward" = "skyblue", "naïf" = "grey"))


#seconde courbe
c_for2 <- rep(0, 8)
c_back2 <- rep(0, 8)

for(i in 2:9){
  c_for2[i - 1] <- median(microbenchmark(forward_discrete(rep(1, 100 * i), mu, p, a)[[1]])$time)
  c_back2[i - 1] <- median(microbenchmark(backward_discrete(rep(1, 100 * i), mu, p, a)[[1]])$time)
}

df2 <- data.frame(x = abscisse*100, y1 = c_for2*10**-6, y2 = c_back2*10**-6)

ggplot(df2, aes(x)) + geom_line(aes(y = y1, colour = "forward"), linewidth = 2) + geom_line(aes(y = y2, colour = "backward"), linewidth = 2) + scale_x_continuous(name = "Longueur de la chaîne") + scale_y_continuous(name = "Temps d'éxectution (ms)") + theme_bw() +
  scale_color_manual(name="", values = c("forward" = "brown1", "backward" = "skyblue"))

#-----------------------------------------------------

#Problematique 2: Test de EM

#Test cas discret Exemple 1

em_discrete(hmm_discrete_8$observed, mu, p, a)
em_discrete(hmm_discrete_1000$observed, mu, p, a)



#Resultat EM vs Gibbs: convergence de EM continu


#Pour sigma = 0.5
sigma_gibbs <- .5

test_b <- continuous_hmm_generator(mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs, 1000)
y_em <- test_b$observed

#Histogramme 1 figure 2
ggplot(data = test_b, aes(x = observed)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, color = "black", fill = "cornflowerblue") +
  labs(x = "Valeur", y = "Densité") + geom_vline(aes(xintercept = mean(observed)), color="black", linetype="dashed", size=1) + theme_bw()

#Initialisation des parametres (Voir EM vs Gibbs)
n <- nrow(p_gibbs)
t_em <- length(y_em)
r_em <- max(y_em) - min(y_em)
mu_em <- rep(1/n, n)
m_em <- min(y_em) + (1:n - 1)*r_em/n + r_em/(2*n)

x_em <- rep(0, t_em)
for(k in 1:t_em){
  x_em[k] <- which.min((y_em[k] - m_em)**2) - 1
}

sigma2_em <- mean((y_em - m_em[x_em + 1])**2)

p_em <- matrix(0, nrow=n, ncol=n)
for (k in 1:n) { 
  n_k <- sapply(0:(n-1), count_transitions, vec = x_em, i = k - 1) 
  p_em[k,] <- n_k / sum(n_k)
}

#Courbe 1 Figure 3
res_em_b <- em_continuous(y_em, mu_em, p_em, m_em, sigma2_em)

m_em_b_1 <- sapply(res_em_b$m, "[[", 1)
m_em_b_2 <- sapply(res_em_b$m, "[[", 2)
m_em_b_3 <- sapply(res_em_b$m, "[[", 3)


ggplot(data = data.frame(m_em_b_1, m_em_b_2, m_em_b_3), aes(x = 1:101)) +
  geom_line(aes(y = m_em_b_1), color = "darkblue", linewidth = 1) +
  geom_line(aes(y = m_em_b_2), color = "darkgreen", linewidth = 1) +
  geom_line(aes(y = m_em_b_3), color = "darkred", linewidth = 1) +
  labs(x = "Iteration de l'algorithme EM", y = "Valeur") +
  theme(plot.title = element_text(size = 4, hjust = .5)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  geom_hline(yintercept = -2, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw()


#Pour sigma = 1
sigma_gibbs <- 1

test_a <- continuous_hmm_generator(mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs, 1000)
y_em <- test_a$observed

#Histogramme 2 figure 2
ggplot(data = test_a, aes(x = observed)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, color = "black", fill = "cornflowerblue") +
  labs(x = "Valeur", y = "Densité") + geom_vline(aes(xintercept = mean(observed)), color="black", linetype="dashed", size=1) + theme_bw()


#Initialisation des parametres (Voir EM vs Gibbs)
n <- nrow(p_gibbs)
t_em <- length(y_em)
r_em <- max(y_em) - min(y_em)
mu_em <- rep(1/n, n)
m_em <- min(y_em) + (1:n - 1)*r_em/n + r_em/(2*n)

x_em <- rep(0, t_em)
for(k in 1:t_em){
  x_em[k] <- which.min((y_em[k] - m_em)**2) - 1
}

sigma2_em <- mean((y_em - m_em[x_em + 1])**2)

p_em <- matrix(0, nrow=n, ncol=n)
for (k in 1:n) { 
  n_k <- sapply(0:(n-1), count_transitions, vec = x_em, i = k - 1) 
  p_em[k,] <- n_k / sum(n_k)
}

#Courbe 2 Figure 3
res_em_a <- em_continuous(y_em, mu_em, p_em, m_em, sigma2_em)

m_em_a_1 <- sapply(res_em_a$m, "[[", 1)
m_em_a_2 <- sapply(res_em_a$m, "[[", 2)
m_em_a_3 <- sapply(res_em_a$m, "[[", 3)


ggplot(data = data.frame(m_em_a_1, m_em_a_2, m_em_a_3), aes(x = 1:101)) +
  geom_line(aes(y = m_em_a_1), color = "darkblue", linewidth = 1) +
  geom_line(aes(y = m_em_a_2), color = "darkgreen", linewidth = 1) +
  geom_line(aes(y = m_em_a_3), color = "darkred", linewidth = 1) +
  labs(x = "Iteration de l'algorithme EM", y = "Valeur") +
  theme(plot.title = element_text(size = 4, hjust = .5)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  geom_hline(yintercept = -2, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw()

#Pour sigma = 1.5
sigma_gibbs <- 1.5

test_c <- continuous_hmm_generator(mu_gibbs, p_gibbs, m_gibbs, sigma_gibbs, 1000)
y_em <- test_c$observed

#Histogramme 3 figure 2
ggplot(data = test_c, aes(x = observed)) +
  geom_histogram(aes(y = ..density..), binwidth = .2, color = "black", fill = "cornflowerblue") +
  labs(x = "Valeur", y = "Densité") + geom_vline(aes(xintercept = mean(observed)), color="black", linetype="dashed", size=1) + theme_bw()

#Initialisation des parametres (Voir EM vs Gibbs)
n <- nrow(p_gibbs)
t_em <- length(y_em)
r_em <- max(y_em) - min(y_em)
mu_em <- rep(1/n, n)
m_em <- min(y_em) + (1:n - 1)*r_em/n + r_em/(2*n)

x_em <- rep(0, t_em)
for(k in 1:t_em){
  x_em[k] <- which.min((y_em[k] - m_em)**2) - 1
}

sigma2_em <- mean((y_em - m_em[x_em + 1])**2)

p_em <- matrix(0, nrow=n, ncol=n)
for (k in 1:n) { 
  n_k <- sapply(0:(n-1), count_transitions, vec = x_em, i = k - 1) 
  p_em[k,] <- n_k / sum(n_k)
}

#Courbe 3 Figure 3
res_em_c <- em_continuous(y_em, mu_em, p_em, m_em, sigma2_em)

m_em_c_1 <- sapply(res_em_c$m, "[[", 1)
m_em_c_2 <- sapply(res_em_c$m, "[[", 2)
m_em_c_3 <- sapply(res_em_c$m, "[[", 3)

ggplot(data = data.frame(m_em_c_1, m_em_c_2, m_em_c_3), aes(x = 1:101)) +
  geom_line(aes(y = m_em_c_1), color = "darkblue", linewidth = 1) +
  geom_line(aes(y = m_em_c_2), color = "darkgreen", linewidth = 1) +
  geom_line(aes(y = m_em_c_3), color = "darkred", linewidth = 1) +
  labs(x = "Iteration de l'algorithme EM", y = "Valeur") +
  theme(plot.title = element_text(size = 4, hjust = .5)) +
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  geom_hline(yintercept = -2, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw()

#-----------------------------------------------------

#Problematique 3: Test Gibbs

#On reprend les chaines generees pour EM

#Pour sigma = 0.5
#Graphe 1 figure 4

res_gibbs_b <- gibbs_sampler(test_b$observed, 3)

m_gibbs_b_1 <- sapply(res_gibbs_b$m, "[[", 1)
m_gibbs_b_2 <- sapply(res_gibbs_b$m, "[[", 2)
m_gibbs_b_3 <- sapply(res_gibbs_b$m, "[[", 3)
m_gibbs_b_df <- data.frame(iter = 1:length(m_gibbs_b_1),
                           m_1 = m_gibbs_b_1,
                           m_2 = m_gibbs_b_2,
                           m_3 = m_gibbs_b_3)

ggplot(data = m_gibbs_b_df, aes(x = iter)) +
  geom_line(aes(y = m_1), color = "darkblue") +
  geom_line(aes(y = m_2), color = "darkgreen") +
  geom_line(aes(y = m_3), color = "darkred") +
  labs(x = "Iteration de Gibbs", y = "Valeur") +
  theme(plot.title = element_text(size = 6)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  geom_hline(yintercept = -2, color = "darkgreen", linetype = "dashed", linewidth = 1) + geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw()

#Pour sigma = 1
#Graphe 2 figure 4

res_gibbs_a <- gibbs_sampler(test_a$observed, 3)


m_gibbs_a_1 <- sapply(res_gibbs_a$m, "[[", 1)
m_gibbs_a_2 <- sapply(res_gibbs_a$m, "[[", 2)
m_gibbs_a_3 <- sapply(res_gibbs_a$m, "[[", 3)
m_gibbs_a_df <- data.frame(iter = 1:length(m_gibbs_a_1),
                           m_1 = m_gibbs_a_1,
                           m_2 = m_gibbs_a_2,
                           m_3 = m_gibbs_a_3)

ggplot(data = m_gibbs_a_df, aes(x = iter)) +
  geom_line(aes(y = m_1), color = "darkblue") +
  geom_line(aes(y = m_2), color = "darkgreen") +
  geom_line(aes(y = m_3), color = "darkred") +
  labs(x = "Iteration de Gibbs", y = "Valeur") +
  theme(plot.title = element_text(size = 6)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  geom_hline(yintercept = -2, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw()

#Pour sigma = 1.5
#Graphe 3 figure 4

res_gibbs_c <- gibbs_sampler(test_c$observed, 3)


m_gibbs_c_1 <- sapply(res_gibbs_c$m, "[[", 1)
m_gibbs_c_2 <- sapply(res_gibbs_c$m, "[[", 2)
m_gibbs_c_3 <- sapply(res_gibbs_c$m, "[[", 3)
m_gibbs_c_df <- data.frame(iter = 1:length(m_gibbs_c_1), 
                           m_1 = m_gibbs_c_1, 
                           m_2 = m_gibbs_c_2, 
                           m_3 = m_gibbs_c_3)

ggplot(data = m_gibbs_c_df, aes(x = iter)) +
  geom_line(aes(y = m_1), color = "darkblue") +
  geom_line(aes(y = m_2), color = "darkgreen") +
  geom_line(aes(y = m_3), color = "darkred") +
  labs(x = "Iteration de Gibbs", y = "Valeur") +
  theme(plot.title = element_text(size = 6)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  geom_hline(yintercept = -2, color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkblue", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw()


#Figure 5 autocorrelation

compute_autocorrelation <- function(x, lag.max) {
  pos_acf <- acf(x, lag.max = lag.max, plot = FALSE)$acf[-1]
  neg_acf <- acf(rev(x), lag.max = lag.max, plot = FALSE)$acf[-1]
  return(c(rev(neg_acf), 1, pos_acf))
}


# Calcul des autocorrelations empiriques dans les estimations de m_2 pour sigma = 1 et 0.5
lag.range <- 30
m_a2_acf <- compute_autocorrelation(m_gibbs_a_2, lag.range)
m_b2_acf <- compute_autocorrelation(m_gibbs_b_2, lag.range)

acf_df <- data.frame(lag = (-lag.range):lag.range, m_a_acf = m_a2_acf, m_b_acf = m_b2_acf)

ggplot(acf_df, aes(lag)) + geom_line(aes(y = m_a_acf), linewidth = 1.1, color = "brown1", linetype = "longdash") + geom_line(aes(y = m_b_acf), linewidth = 1.1, color = "skyblue") + scale_x_continuous(name = "Sweep lag") + scale_y_continuous(name = "Autocorrelation") + theme_bw()

#-----------------------------------------------------

#Problematique 3: Test Viterbi

#Test sur Exemple 1

#chaine de taille 8
y <- c(0, 1, 1, 1, 1, 0, 1, 0)

x_opti <- viterbi_discrete(y, mu, p, a)$path
proba_opti <- viterbi_discrete(y, mu, p, a)$prob

#Figure 6 graphe 1
chain_possible <- as.matrix(expand.grid(rep(list(states),length(y))))

result <- rep(0, 256)

for (i in 1 : 256){
  result[i] <- proba(y, chain_possible[i,], mu, p, a)
}

data_result <- data.frame(probs = result, number = 1:256)

ggplot(data = data_result, aes(x = number, y=probs)) +
  geom_bar(stat = "identity", fill = "cornflowerblue")+
  theme_bw() + labs(x = "", y = "Probabilité jointe")


#Test sur nouvelle exemple plus complexe

#chaine de taille 6 pour un nouveau modele

mu_viterbi <- c(.4, .3, .3)
a_viterbi <- matrix(c(.37, .45, .18, .26, .70, .04, .33, .45, .22), nrow = 3, ncol = 3, byrow = T)
p_viterbi <- matrix(c(.5, .3, .2, .2, .8, 0, .3, .3, .4), nrow = 3, ncol = 3, byrow = T)
states_viterbi <- c(0, 1, 2)

y2 <- c(1, 1, 0, 1, 1, 1)

x_opti <- viterbi_discrete(y2, mu_viterbi, p_viterbi, a_viterbi)$path
proba_opti <- viterbi_discrete(y2, mu_viterbi, p_viterbi, a_viterbi)$prob

#Figure 6 graphe 2

chain_possible <- as.matrix(expand.grid(rep(list(states_viterbi), length(y2))))

taille <- 3**6
result <- rep(0, taille)
for (i in 1:taille){
  result[i]<-proba(y2, chain_possible[i,], mu_viterbi, p_viterbi, a_viterbi)
}

data_result <- data.frame(probs = result, number = 1:taille)

ggplot(data = data_result, aes(x = number, y = probs)) +
  geom_bar(stat = "identity", fill = "cornflowerblue") +
  theme_bw() + labs(x = "", y = "Probabilité jointe")

#Il faut ouvrir en version zoomee dans une nouvelle fenetre le graphique pour voir la chaine optimale


#Comparaison complexite naif et Viterbi
#non present sur le manuscrit pour eviter une redite avec la partie 1


c_viterbi <- rep(0, 10)
c_naif2 <- rep(0, 10)
for(i in 2 : 11){
  c_viterbi[i - 1] <- median(microbenchmark(viterbi_discrete(rep(1, i), mu, p, a))$time)
  c_naif2[i - 1] <- median(microbenchmark(naif_viterbi_discrete(rep(1, i), mu, p, a, states))$time)
}

#comparaison naif viterbi
abscisse3<- 2 : 11
df3 <- data.frame(x = abscisse3, y1 = c_viterbi * 10**-6, y2 = c_naif2 * 10**-6)
ggplot(df3, aes(x)) + geom_line(aes(y = y1, colour = "Viterbi"), linewidth = 2) + geom_line(aes(y = y2, colour = "Naïf"), linewidth = 2) + scale_x_continuous(name = "Longueur de la chaine") + scale_y_continuous(name = "Temps d'éxectution (ms)") +
  theme_bw() + scale_color_manual(name="", values = c("Viterbi" = "brown1", "Naïf" = "grey"))


#efficacite de viterbi pour de longues chaines
c_viterbi2 <- rep(0, 10)
for(i in 2 : 11){
  c_viterbi2[i - 1] <- median(microbenchmark(viterbi_discrete(rep(1, 100 * i), mu, p, a))$time)
}

df4 <- data.frame(x = abscisse3 * 100, y1 = c_viterbi2 * 10**-6)
ggplot(df4, aes(x)) + geom_line(aes(y = y1, colour = "Viterbi"), linewidth = 2) + scale_x_continuous(name="Longueur de la chaine") + scale_y_continuous(name = "Temps d'éxectution (ms)") + 
  theme_bw() +  scale_color_manual(name="", values = c("Viterbi" = "brown1"))

