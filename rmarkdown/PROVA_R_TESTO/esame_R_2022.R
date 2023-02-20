#---------------------------------------------#
# ESAME DI STATISTICA - ING. MATEMATICA #
# Prova di laboratorio R - 05/07/2022 #
#---------------------------------------------#

#L'abalone e' un mollusco pregiato, dotato di una sola conchiglia e diffuso in varie aree del globo. Sono stati raccolti 200 esemplari attorno ad alcune isole dello stretto di Bass, che separa l'Australia dalla Tasmania, registrandone il sesso (M=Maschio, F=Femmina o I=Infante), la lunghezza della conchiglia (ovvero la misura maggiore) in mm, il diametro in mm e l'altezza in mm. Questi dati sono riportati nel file abalone4.csv
#Rispondete alle seguenti domande, riportando il codice R, le immagini dei plot ove richiesti e commentando i risultati ottenuti
# 0) Importare il dataset.
# 1) Calcolare le frequenze assolute e relative e la moda delle variabili categoriche e descrivere queste variabili con uno o piu' grafici opportuni.
# 2) Calcolare media, mediana, quartili, minimo, massimo, 10° percentile e 90° percentile delle variabili numeriche.
# 3) Verificare graficamente e numericamente la correlazione tra le variabili Length e Height.
# 4) Analizzare gli scatter plot tra tutte le coppie di variabili numeriche, colorando i punti in base al sesso.
# 5) Analizzare la distribuzione della variabile Height con gli strumenti grafici piu' opportuni, sia globale sia rispetto ai diversi valori della variabile Sex.
# 6) Verificare se è possibile affermare che la varianza della variabile Height e' maggiore di 52 a un livello di significativita' del 5%.
# 7) Testare se e' possibile affermare che i maschi (M) e le femmine (F) hanno altezze diverse. Assumere uguali varianze tra i due gruppi.
# 8) Calcolare un intervallo di confidenza bilatero al 99% per la differenza delle medie tra i gruppi M e I (altezza_maschi - altezza_infanti).
# 9) Calcolare la potenza del test del punto 7) via MC, per una differenza delle medie mu.M - mu.F = 2, usando 1000 iterazioni. Generare i campioni normali con la stessa numerosita' n=200 e la stessa deviazione standard = 6 per entrambi i gruppi.
# 10) Mostrare la curva di potenza del test del punto 7) al variare della differenza delle medie tra -10 e +10; usare ancora MC per il calcolo della potenza, in modo analogo al punto 9, calcolandola con un passo di 0.2 per la differenza delle medie.

setwd("C:/Users/vitto/OneDrive - Politecnico di Milano/statistica/MATERIALE 21-22")

rm(list=ls())
graphics.off()

#0) Importare il dataset
dati = read.csv("abalone4.csv")
dati
colnames(dati)
dim(dati)
n=dim(dati)[1]
attach(dati)

# 1) Calcolare le frequenze assolute e relative e la moda delle variabili categoriche e descriverle con un grafico opportuno

#frequenze assolute
freq.assolute = table(Sex)
freq.assolute
#frequenza relative
freq.relative = prop.table(table(Sex))
freq.relative
#la moda e' I
freq.assolute[freq.assolute==max(freq.assolute)]

dev.new()
pie(freq.assolute, col=c('cyan', 'coral', 'chartreuse'), main = 'Sex')

dev.new()
barplot(freq.assolute, col=c('cyan', 'coral', 'chartreuse'), main = 'Sex')


# 2) Calcolare media, mediana, quartili, minimo, massimo, 10° percentile e 90° percentile delle variabili numeriche
summary(Length)
quantile(Length, c(0.10, 0.90), type=2)


summary(Height)
quantile(Height, c(0.10, 0.90), type=2)


summary(Diameter)
quantile(Diameter, c(0.10, 0.90), type=2)

# 3) Verificare graficamente e numericamente la correlazione tra le variabili Length e Height

dev.new()
plot(Length, Height, pch=20, col='blue')

cor(Length, Height)

#COMMENTO: sia graficamente che numericamente si evince un'elevata correlazione

# 4) Analizzare gli scatter plot tra tutte le coppie di variabili numeriche, colorando i punti in base al sesso
colors = rep('cyan', n)
colors[Sex=='I'] = 'coral'
colors[Sex=='M'] = 'chartreuse'
dati.numerici = dati[,-1]
dev.new()
pairs(dati.numerici, col=colors, pch=20)

#COMMENTO: tra Length e Diameter si evince una correlazione ancora più alta di quella tra Length e Height, e anche tra Diameter e Height sembra esserci una buona correlazione. Gli Infanti (I) tendono ad avere valori più bassi, mentre non risultano evidenti differenze tra M e F.   

# 5) Analizzare la distribuzione della variabile Height con gli strumenti grafici piu' opportuni, sia globale sia rispetto ai diversi valori della variabile Sex

dev.new()
hist(Height)
#Distribuzione unimodale, con una leggera asimmetria negativa

dev.new()
par(mfrow=c(3,1))
hist(Height[Sex=='I'], col='cyan')
hist(Height[Sex=='F'], col='coral')
hist(Height[Sex=='M'], col='chartreuse')

dev.new()
boxplot(Height~Sex, col=c('cyan', 'coral', 'chartreuse'))

#La mediana nei gruppi M e F e' simile, la varianza e' leggermente piu' alta nel gruppo M, che ha una distribuzione meno simmetrica, con una coda superiore piu' pesante. Il gruppo I presenta una mediana decisamente piu' bassa rispetto agli altri due gruppi e una varianza leggermente piu' elevata rispetto al gruppo M.

# 6) Verificare se è possibile affermare che la varianza della variabile Height e' maggiore di 52 a un livello di significativita' del 5%

# verifichiamo la normalità dei dati
qqnorm(Height)
qqline(Height,lwd=2,col='red')


shapiro.test(dati$Height)

#il p-value dello shapiro test è alto e anche il qqplot mostra un andamento normale
#procediamo quindi con il test chi-quadro


alpha = 0.05
sigma2_0 = 52
library(EnvStats)
res = varTest(Height, alternative="greater", conf.level = 1-alpha, sigma.squared = sigma2_0)

res$p.value
#al 5% non possiamo affermare che la varianza è maggiore di 52


# 7) Testare se e' possibile affermare che i maschi (M) e le femmine (F) hanno altezze diverse. Assumere uguali varianze tra i due gruppi


height.m = Height[Sex=='M']
height.f = Height[Sex=='F']

#Verifichiamo la normalita' dei due gruppi
qqnorm(height.f)
qqline(height.f, lwd=2, col='red')

qqnorm(height.m)
qqline(height.m, lwd=2, col='red')

shapiro.test(height.f)
shapiro.test(height.m)

#l'assunzione di normalita' e' ragionevole, visti i p-value dello shapiro test e i qqplot

#var.test(height.m, height.f, alternative = "two.sided", ratio = 1)
#(Il test di confronto delle varianze ha pvalue 0.15, quindi l'assunzione di uguaglianza delle varianze fatta nel testo e' ragionevole)

res = t.test(height.m, height.f, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 1 - alpha)
res$p.value

#dato il valore del p-value (0.50), non è possibile affermare che i due gruppi abbiano medie differenti

# 8) Calcolare un intervallo di confidenza bilatero al 99% per la differenza delle medie tra i gruppi M e I (altezza_maschi - altezza_infanti)

#verifichiamo la normalita' del gruppo I (quella del gruppo M e' gia' stata verificata al punto precedente)
height.I = Height[Sex=='I']
qqnorm(height.I)
qqline(height.I, lwd=2, col='red')

shapiro.test(height.I)

# Dato il p-value dello shapiro test, possiamo assumere la normalita'

alpha = 0.01
res = t.test(height.m, height.I, alternative = "two.sided", conf.level = 1 - alpha, mu = 0, var.equal = TRUE, paired = FALSE)
res$conf.int

#(5.638539, 11.447999) e' l'IC bilatero al 99% per la differenza h.M - h.I

# 9) Calcolare la potenza del test del punto 7) via MC, per una differenza delle medie mu.M - mu.F = 2, usando 1000 iterazioni. Generare i campioni normali con la stessa numerosita' n=200 e la stessa deviazione standard = 6

M = 1000
n = 200
dev.std = 6
alpha = 0.05
esito = rep(0, M)
for(i in 1:M){
  m = rnorm(n, mean=2, sd=dev.std)
  f = rnorm(n, mean=0, sd=dev.std)
  res = t.test(m, f, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 1 - alpha)
  if(res$p.value < alpha)
    esito[i] = 1
}
power = mean(esito)
power
#La potenza stimata e' di 0.91, elevata

# 10) Mostrare la curva di potenza del test del punto 7) al variare della differenza delle medie tra -10 e +10; usare ancora MC per il calcolo della potenza, in modo analogo al punto 9), calcolandola con un passo di 0.2 per la differenza delle medie.

valori = seq(-10, 10, 0.2)
n.valori = length(valori)
power = rep(0, n.valori)
for(j in 1:n.valori){
  esito = rep(0, M)
  for(i in 1:M){
    m = rnorm(n, mean=valori[j], sd=dev.std)
    f = rnorm(n, mean=0, sd=dev.std)
    res = t.test(m, f, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 1 - alpha)
    if(res$p.value < alpha)
      esito[i] = 1
  }
  power[j] = mean(esito)
}

dev.new()
plot(valori, power, type='l', col='blue')

#La curva è simmetrica rispetto a zero, come e' ragionevole aspettarsi da un test bilatero, e raggiunge valori elevati di potenza a partire da una differenza di circa 2
