library(Familias)
p = 1
a = "1"
nn = "L1"
model ="Equal"
r = 1
a=FamiliasLocus(frequencies = p, allelenames= a, name = nn,
              MutationModel = model,
              MutationRate  = r)

model = "stepwise"
als= c("16", "16.1", "17", "17.1", "18")
als = 10:14
nall = length(alleles)
p = rep(1/5,5) # not used
names(afreq) = alleles
r = 0.003
r2 = 0.001
Range = 0
M=mutationMatrix(model = "stepwise",
               matrix = NULL, alleles = als, afreq = p,rate = r,
               rate2 = r2, range = Range)
M
library(Familias)
M2 = FamiliasLocus(frequencies = p, allelenames = als, name ="L1",
                   MutationModel = "Stepwise", MutationRate = r,
                   MutationRate2 = r2, MutationRange = Range)$femaleMutationMatrix
all(M1==M2)
validateMutationMatrix(M1)
