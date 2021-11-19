CDFt <- function(ObsRp, DataGp, DataGf, npas = 1000, dev = 2){
	dyn.load("CDFt.so")
	output <- .Call(
		"cdft",
		as.numeric(na.omit(ObsRp)),
		as.numeric(na.omit(DataGp)),
		as.numeric(na.omit(DataGf)),
		as.integer(npas),
		as.numeric(dev)
	)
	names(output) <- c("x", "FRp", "FGp", "FGf", "FRf", "DS")
	return(output)
}

n <- 1000000
vrp <- rnorm(n)
vgp <- sample(rnorm(n) + 2)
vgf <- sample(rnorm(n) + 4)

print("-----------")
print( 
       system.time({
	       a <- CDFt(vrp, vgp, vgf, npas = 20)
       })
)
print("-----------")
print( 
       system.time({
	       b <- CDFt::CDFt(vrp, vgp, vgf, npas = 20)
       })
)
print(all.equal(a, b))
print(all.equal(a$DS, b$DS))
# plot(a$FGp, b$FGp - a$FGp)
# plot(a$FGf, b$FGf - a$FGp)
# plot(a$DS, b$DS)

# Rprof("boot.out")
# a <- CDFt(vrp, vgp, vgf, npas = 20)
# Rprof(NULL)
# Rprof("boot2.out")
# b <- CDFt::CDFt(vrp, vgp, vgf, npas = 20)
# Rprof(NULL)
