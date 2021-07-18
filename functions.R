make_corr <- function(D, 
                      phi, 
                      form = "exponential") {
    # evaluate the kernel/basis in 1d
    if (!(form %in% c("exponential", "gaussian"))) {
        stop('form must be either "exponential" or "gaussian" ...')
    }
    
    if (form == "exponential") { 
        corr_fun <- exp(-D / phi)
    } 
    if (form == "gaussian") { 
        corr_fun <- exp(-D^2 / phi)
    }
    
    return(corr_fun)
        
}

mcmc_gp <- function(y, 
                    X, 
                    D,
                    form = "exponential",
                    tuning_parameters,
                    n_mcmc = 5000, 
                    burnin = 2500, 
                    n_message = 500) {
    library(mvnfast)
    library(truncnorm)
    
    y <- c(y)
    N <- length(y)
    p <- ncol(X)
    
    ## get prior parameters
    a_sigma <- 0.01
    b_sigma <- 0.01
    a_tau <- 0.01
    b_tau <- 0.01
    Sigma_beta <- 100 * diag(p)
    Sigma_beta_inv <- 1/100 * diag(p)
    
    ## initialize parameters
    sigma2<- runif(1, 1, 5)
    tau2<- runif(1, 1, 5)
    beta <- matrix(rnorm(p, 0, 1), nrow = p, ncol = 1)
    phi <- runif(1, 0, 1)
    C <- make_corr(D, phi, form = "exponential")
    
    ## precalculate values
    Xbeta <- c(X %*% beta)
    I_p <- diag(p)
    I_N <- diag(N)
    
    ## set up save variables
    sigma2_save <- rep(0, n_mcmc-burnin)
    tau2_save <- rep(0, n_mcmc-burnin)
    phi_save <- rep(0, n_mcmc-burnin)
    beta_save <- matrix(0, n_mcmc-burnin, p)
    
    ## tuning and acceptance variables
    phi_tune <- tuning_parameters$phi_tune
    phi_accept <- 0
    phi_accept_batch <- 0
    
    sigma2_tune <- tuning_parameters$sigma2_tune
    sigma2_accept <- 0
    sigma2_accept_batch <- 0
    
    tau2_tune <- tuning_parameters$tau2_tune
    tau2_accept <- 0
    tau2_accept_batch <- 0
    
    ## mcmc loop
    for(k in 1:n_mcmc) {
        # if (k%% n_message == 0) {
        #     message("Iteration ", k , " out of ", n_mcmc)
        # }
        
        Sigma_y <- sigma2 * diag(N) + tau2*C
        Sigma_y_inv <- solve(Sigma_y)
        Sigma_y_chol <- chol(Sigma_y)
        
        ## sample beta ---------------------------------------------------------
        A_inv <- solve(t(X) %*% Sigma_y_inv %*% X + Sigma_beta_inv)
        b <- t(X) %*% Sigma_y_inv %*% y
        beta <- c(rmvn(1, A_inv %*% b, A_inv))
        Xbeta <- c(X %*% beta)
        
        ## sample phi ----------------------------------------------------------
        phi_star <- rtruncnorm(1, a=0, b=Inf, phi, phi_tune)
        C_star <- make_corr(D, phi_star, form = "exponential")
        Sigma_y_star <- sigma2 * I_N + tau2 * C_star
        Sigma_y_chol_star <- chol(Sigma_y_star)
        
        mh1 <- mvnfast::dmvn(y, Xbeta, Sigma_y_chol_star, isChol=TRUE, log=TRUE) + 
            log(dtruncnorm(phi, a=0, b=Inf, phi_star, phi_tune))
        
        mh2 <- mvnfast::dmvn(y, Xbeta, Sigma_y_chol, isChol=TRUE, log=TRUE) + 
            log(dtruncnorm(phi_star, a=0, b=Inf, phi, phi_tune))
        
        mh <- exp(mh1 - mh2) 
        if (mh > runif(1)) {
            phi <- phi_star
            C <- C_star
            Sigma_y <- Sigma_y_star
            Sigma_y_chol <- Sigma_y_chol_star
            phi_accept <- phi_accept + 1 / n_mcmc
            phi_accept_batch <- phi_accept_batch + 1 / 50
        }
        
        #### automatically update tuning for phi
        if (k %% 50 == 0) {
            delta = 1.0 / sqrt(k)
            if (phi_accept_batch > 0.44) {
                phi_tune <- exp(log(phi_tune) + delta)
            } else {
                phi_tune <- exp(log(phi_tune) - delta)
            }
            phi_accept_batch <- 0
        }
        
        ##sample sigma2 --------------------------------------------------------
        sigma2_star <- rtruncnorm(1, a=0, b=Inf, mean = sigma2, sd = sigma2_tune)
        Sigma_y_star <- sigma2_star * I_N + tau2 * C
        Sigma_y_chol_star <- chol(Sigma_y_star)
        
        mh1 <- mvnfast::dmvn(y, Xbeta, Sigma_y_chol_star, isChol=TRUE, log=TRUE) + 
            (dgamma(sigma2_star, a_sigma, b_sigma,
                   log=TRUE)) + 
            log(dtruncnorm(sigma2, a=0, b=Inf, mean = sigma2_star,
                           sd = sigma2_tune))
        
        mh2 <- mvnfast::dmvn(y, Xbeta, Sigma_y_chol, isChol=TRUE, log=TRUE) + 
            (dgamma(sigma2, a_sigma, a_sigma,
                   log=TRUE)) + 
            log(dtruncnorm(sigma2_star, a=0, b=Inf, sigma2, 
                           sigma2_tune))
        
        mh <- exp(mh1 - mh2) 
        if (mh > runif(1)) {
            sigma2 <- sigma2_star
            Sigma_y <- Sigma_y_star
            Sigma_y_chol <- Sigma_y_chol_star
            sigma2_accept <- sigma2_accept + 1 / n_mcmc
            sigma2_accept_batch <- sigma2_accept_batch + 1 / 50
        }
        
        ## automatically update tuning for sigma2
        if (k %% 50 == 0) {
            delta = 1.0 / sqrt(k)
            if (sigma2_accept_batch > 0.44) {
                sigma2_tune <- exp(log(sigma2_tune) + delta)
            } else {
                sigma2_tune <- exp(log(sigma2_tune) - delta)
            }
            sigma2_accept_batch <- 0
        }
        
        ## sample tau2 ---------------------------------------------------------
        tau2_star <- rtruncnorm(1, a=0, b=Inf, mean = tau2, sd = tau2_tune)
        Sigma_y_star <- sigma2 * I_N + tau2_star * C
        Sigma_y_chol_star <- chol(Sigma_y_star)
        
        mh1 <- mvnfast::dmvn(y, Xbeta, Sigma_y_chol_star, isChol=TRUE, log=TRUE) + 
            (dgamma(tau2_star, a_tau, b_tau, log=TRUE)) + 
            log(dtruncnorm(tau2, a=0, b=Inf, tau2_star, tau2_tune))
        
        mh2 <- mvnfast::dmvn(y, Xbeta, Sigma_y_chol, isChol=TRUE, log=TRUE) + 
            (dgamma(tau2, a_tau, a_tau, log=TRUE)) + 
            log(dtruncnorm(tau2_star, a=0, b=Inf, tau2, tau2_tune))
        
        mh <- exp(mh1 - mh2) 
        if (mh > runif(1)) {
            tau2 <- tau2_star
            Sigma_y <- Sigma_y_star
            Sigma_y_chol <- Sigma_y_chol_star
            tau2_accept <- tau2_accept + 1 / n_mcmc
            tau2_accept_batch <- tau2_accept_batch + 1 / 50
        }
        
        ## automatically update tuning for sigma2_epsilon
        if (k %% 50 == 0) {
            delta = 1.0 / sqrt(k)
            if (tau2_accept_batch > 0.44) {
                tau2_tune <- exp(log(tau2_tune) + delta)
            } else {
                tau2_tune <- exp(log(tau2_tune) - delta)
            }
            tau2_accept_batch <- 0
        }
        
        ## save MCMC variables
        if(k > burnin) {
            i <- k-burnin
            beta_save[i, ] <- beta
            phi_save[i] <- phi
            sigma2_save[i] <- sigma2
            tau2_save[i] <- tau2
        }
    }
    
    return(
        list(
            beta = beta_save,
            phi = phi_save,
            sigma2 = sigma2_save,
            tau2 = tau2_save,
            phi_accept = phi_accept,
            sigma2_accept = sigma2_accept,
            tau2_accept = tau2_accept,
            n_samples = n_mcmc-burnin
        )
    )
}

sample_eta <- function(y,
                       mcmc_output,
                       D
                       ) {
    N <- dim(D)[1]
    n_samples <- mcmc_output$n_samples
    eta <- matrix(0, n_samples, N) 
    I_N <- diag(N)
    
    for (k in 1:n_samples) {
        # if (k %% 500 == 0) {
        #     message("Iteration ", k , " out of ", n_samples)
        # }
        ## sample z
        sigma2_inv <- 1/mcmc_output$sigma2[k] * I_N
        C <- make_corr(D, mcmc_output$phi[k], form = "exponential")
        
        A_inv <- solve(sigma2_inv + solve(C / mcmc_output$tau2[k]))
        b <- y / mcmc_output$sigma2[k]
        eta[k, ] <- c(mvnfast::rmvn(1, A_inv %*% b, A_inv)) 
    }
    
    mcmc_output$eta <- eta
    
    return(mcmc_output)
}