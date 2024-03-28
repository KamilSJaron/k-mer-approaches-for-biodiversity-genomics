library('GenomeTelescope')

##### BEGONIAS (A)

begonia <- read.table('Bjohnst_K21_pre_filtering.histo.txt', header = F, col.names = c('cov', 'freq'))
# begonia <- read.table('Bjohnst_k21_post_filtering.histo', header = F, col.names = c('cov', 'freq'))

head(begonia)

# coverage_barplot(begonia$freq[cov_range], begonia$cov[cov_range])
# lines(3.870e+01)
# 1.557e+02

cov_range <- 1:300
font_size = 1
width = 0.5
bar_heights <- begonia$freq[1:300]
bar_positions <- begonia$cov[1:300]

pdf('begonia_composite_plot.pdf')
plot(bar_heights, type = "n", xlab = "Coverage", ylab = "Frequency", 
    ylim = c(0, 5e6), xlim = range(0:300), 
    cex.lab = font_size, cex.axis = font_size, cex.main = font_size, 
    cex.sub = font_size)
for (i in 1:length(bar_heights)) {
    rect(bar_positions[i] - width, 0, bar_positions[i] + 
        width, bar_heights[i], col = "deepskyblue", border = F)
}

cov_range_fit <- 20:300
x <- begonia$cov[cov_range_fit]
y <- begonia$freq[cov_range_fit]
lengthEst <- 1e9
hetEst <- 0.8
kmerEst <- 50
kmerEst2 <- 120
# TODO: adjust to fit the three peaks all together
# (x, y, kmerEst, lengthEst, hetEst = 0.6)
begonia_model <- nlsLM(y ~ ((het * dnbinom(x, size = kmercov/bias , mu = kmercov )) + 
                            (het * dnbinom(x, size = kmercov2/bias, mu = kmercov2)) + 
                            ((1 - het) * dnbinom(x, size = (kmercov + kmercov2)/bias, mu = kmercov + kmercov2))) * length,
                        start = list(kmercov = kmerEst, kmercov2 = kmerEst2, 
                        bias = 0.5, length = lengthEst, het = hetEst), 
                        control = list(minFactor = 1e-12, maxiter = 40))

# lines(predict(begonia_model) ~ cov_range_fit, lwd = 3, col = 1)

model_env <- begonia_model$m$getEnv()

left_peak <- model_env$het * dnbinom(model_env$x, size = model_env$kmercov/model_env$bias , mu = model_env$kmercov) * model_env$length
right_peak <- model_env$het * dnbinom(model_env$x, size = model_env$kmercov2/model_env$bias , mu = model_env$kmercov2 ) * model_env$length
joint_peak <- (1 - model_env$het) * dnbinom(model_env$x, size = (model_env$kmercov + model_env$kmercov2)/model_env$bias, mu = model_env$kmercov + model_env$kmercov2) * model_env$length

hybrid_col <- "chocolate"
target_col <- "darkorchid4"

lines(left_peak ~ cov_range_fit, lwd = 3, col = hybrid_col)
lines(right_peak ~ cov_range_fit, lwd = 3, col = target_col)
lines(joint_peak ~ cov_range_fit, lwd = 3, col = 1)

legend("topright", c("kmer histogram", "Begonia sp.", 
    "Begonia johnstonii", "Both Begonias"), col = c("deepskyblue", hybrid_col, 
    target_col, "black"), lty = c(NA, 1, 1, 1), lwd = 3, pch = c(15, 
    NA, NA, NA), bty = "n")
dev.off()

##### CORAL (B)
require('GenomeTelescope')

# the spectra can be fetched from
# https://tolqc.cog.sanger.ac.uk/asg/jellyfish/Pocillopora_grandis/genomic_data/jaPocGran1/illumina/kmer/k31/jaPocGran1.k31.hist.txt
coral <- read.table('jaPocGran1.k31.hist.txt', header = F, col.names = c('cov', 'freq'))

cov_range <- 4:400
coverage_barplot(coral$freq[cov_range], coral$cov[cov_range])

cov_range <- 1:400
font_size = 1
width = 0.5
bar_heights <- coral$freq[cov_range]
bar_positions <- coral$cov[cov_range]

pdf('coral_composite_plot.pdf')
plot(bar_heights, type = "n", xlab = "Coverage", ylab = "Frequency", 
    ylim = c(0, 35e6), xlim = range(cov_range), 
    cex.lab = font_size, cex.axis = font_size, cex.main = font_size, 
    cex.sub = font_size)
for (i in 1:length(bar_heights)) {
    rect(bar_positions[i] - width, 0, bar_positions[i] + 
        width, bar_heights[i], col = "deepskyblue", border = F)
}

nlsLM_2peak_model <- function(x, y, estKmercov, estLength, estR, k = 31){
	nlsLM(y ~ (((2*(1-(1-r)^k)) * dnbinom(x, size = kmercov   / bias, mu = kmercov)) +
			   ((1-r)^k)         * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)) * length,
		  start = list(kmercov = estKmercov, bias = 0.5, length = estLength, r = estR),
		  control = list(minFactor=1e-12, maxiter=40))
}

cov_range_1 <- 6:80
k <- 31
cobiont_1 <- nlsLM_2peak_model(coral$cov[cov_range_1], coral$freq[cov_range_1], 15, 100e6, 0.05)

x <- cobiont_1$m$getEnv()$x
y <- cobiont_1$m$getEnv()$y
cov_1n <- coef(cobiont_1)["kmercov"]
het <- coef(cobiont_1)["kmercov2"]
genomesize <- coef(cobiont_1)["length"]
cobiont1_col <- "darkorchid4"

lines(predict(cobiont_1) ~ x, col = cobiont1_col, lwd = 3)


cov_range_2 <- 90:400
cobiont_2 <- nlsLM_2peak_model(coral$cov[cov_range_2], coral$freq[cov_range_2], 150, 100e6, 0.04)

cobiont2_col <- "chocolate"
lines(predict(cobiont_2) ~ cov_range_2, col = cobiont2_col, lwd = 3)

legend("topright", c("kmer histogram", "cobiont", 
    "stone coral"), col = c("deepskyblue", cobiont1_col, 
    cobiont2_col), lty = c(NA, 1, 1), lwd = 3, pch = c(15, 
    NA, NA), bty = "n")
dev.off()