run_example <- function() {
    library(ggplot2)
    source("functions.R")
    set.seed(23)
    
    theme_set(theme_minimal())
    
    ## Simulate data ---------------------------------------------------------------
    N <- 200
    s <- seq(0, 1, length.out = N)
    D <- fields::rdist(s)
    p <- 2 # number of regression covariates 
    
    X <- matrix(mvnfast::rmvn(p, mu = rep(0, N), sigma = diag(N)), nrow = N, ncol = p)
    beta <- matrix(rnorm(p, 0, 1), nrow = p, ncol = 1)
    X_beta <- X %*% beta
    
    phi <- 1/8
    tau2 <- 3^2
    C <- make_corr(D, phi, form = "exponential")
    eta <- t(mvnfast::rmvn(1, mu=rep(0, N), sigma = tau2*C))
    
    sigma2 <- 0.5^2
    
    y <- X_beta + eta + rnorm(N, 0, sqrt(sigma2))
    
    dat <- data.frame(
        s = s,
        z = X_beta + eta,
        y = y
    )
    
    ggplot(data = dat) + 
        geom_point(aes(x=s, y=y)) + 
        geom_line(aes(x=s, y=z))
    
    tuning_parameters <- list(
        phi_tune = 0.2,
        sigma2_tune = 0.03,
        tau2_tune = 2
    )
    
    fit <- mcmc_gp(y, X, D,
                   form = "exponential",
                   tuning_parameters,
                   n_mcmc = 5000, burnin = 2500, n_message = 500)
    
    layout(matrix(1:4, 2, 2)) 
    plot(fit$sigma2, type='l',
         main=paste0("Accept = ", round(fit$sigma2_accept, digits=2)))
    abline(h=sigma2, col="red") 
    plot(fit$tau2, type='l',
         main=paste0("Accept = ", round(fit$tau2_accept, digits=2)))
    abline(h=tau2, col="red") 
    plot(fit$phi, type='l',
         main=paste0("Accept = ", round(fit$phi_accept,digits=2))) 
    abline(h=phi, col="red")
    
    fit <- sample_eta(y, fit, D)
    
    df <- data.frame(x = s,
                     y = y)
    
    df_pred <- data.frame(
        x = s,
        z = X_beta + eta,
        lower_50=apply(t(X %*% t(fit$beta)) + fit$eta, 2, quantile, prob=0.25),
        upper_50=apply(t(X %*% t(fit$beta)) + fit$eta, 2, quantile, prob=0.75), 
        lower_95=apply(t(X %*% t(fit$beta)) + fit$eta, 2, quantile, prob=0.025), 
        upper_95=apply(t(X %*% t(fit$beta)) + fit$eta, 2, quantile, prob=0.975))
    
    ggplot(data=df, aes(x=x, y=y)) + 
        geom_point() +
        geom_line(data=df_pred, aes(x=x, y=z), col="red") + 
        geom_ribbon(data=df_pred, aes(x=x, ymin=lower_50, ymax=upper_50),
                    fill="blue", alpha = 0.5, inherit.aes = FALSE) + 
        geom_ribbon(data=df_pred, aes(x=x, ymin=lower_95, ymax=upper_95),
                    fill="blue", alpha = 0.25, inherit.aes = FALSE)
}
