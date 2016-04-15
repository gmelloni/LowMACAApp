cooccur <- function (mat, type = "spp_site", thresh = TRUE, spp_names = FALSE, 
    true_rand_classifier = 0.1, prob = c("hyper" , "comb"), only_effects = FALSE, 
    eff_standard = TRUE, eff_matrix = FALSE) 
{
    prob <- prob[1]
    if (type == "spp_site") {
        spp_site_mat <- mat
    }
    if (type == "site_spp") {
        spp_site_mat <- t(mat)
    }
    if (spp_names == TRUE) {
        spp_key <- data.frame(num = 1:nrow(spp_site_mat), spp = row.names(spp_site_mat))
    }
    spp_site_mat[spp_site_mat > 0] <- 1
    nsite <- ncol(spp_site_mat)
    nspp <- nrow(spp_site_mat)
    spp_pairs <- choose(nspp, 2)
    obs_cooccur <- prob_cooccur <- exp_cooccur <- matrix(nrow = spp_pairs, 
        ncol = 3)
    #########incidence <- rowSums(spp_site_mat, na.rm = T)
    #########prob_occur <- incidence/nsite
    prob_occur <- incidence <- cbind(1:nspp, rowSums(spp_site_mat, 
        na.rm = T))
    prob_occur[ , 2] <- prob_occur[ , 2]/nsite
    # prob_occur <- cbind(1:nspp, rowSums(spp_site_mat, 
    #     na.rm = T)/nsite)
    #pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), 
    #    style = 3)

    spp_pairs_mat <- t(combn(nspp,2))
    pairs_observations <- apply(spp_pairs_mat, 1, function(x) {
        spp <- x[1]
        spp_next <- x[2]
        pairs <- sum( spp_site_mat[spp, ] * spp_site_mat[spp_next, ] )
        pairs_prob <- sum( prob_occur[spp, 2] * prob_occur[spp_next, 2] )
        pairs_exp <- pairs_prob * nsite
        return(c(pairs, pairs_prob, pairs_exp))
        })

    obs_cooccur <- cbind(spp_pairs_mat, pairs_observations[1,])
    prob_cooccur <- cbind(spp_pairs_mat, pairs_observations[2,])
    exp_cooccur <- cbind(spp_pairs_mat, pairs_observations[3,])

    if (thresh == TRUE) {
        n_pairs <- nrow(prob_cooccur)
        prob_cooccur <- prob_cooccur[exp_cooccur[, 3] >= 1, ,drop=FALSE]
        obs_cooccur <- obs_cooccur[exp_cooccur[, 3] >= 1, ,drop=FALSE]
        exp_cooccur <- exp_cooccur[exp_cooccur[, 3] >= 1, ,drop=FALSE]
        n_omitted <- n_pairs - nrow(prob_cooccur)
        #pb <- txtProgressBar(min = 0, max = (nspp + nrow(obs_cooccur)), 
        #    style = 3)
    }
    output <- t(sapply(1:nrow(obs_cooccur) , function(row) {
        sp1 <- obs_cooccur[row, 1]
        sp2 <- obs_cooccur[row, 2]
        sp1_inc <- incidence[incidence[, 1] == sp1, 2]
        sp2_inc <- incidence[incidence[, 1] == sp2, 2]
        max_inc <- max(sp1_inc, sp2_inc)
        min_inc <- min(sp1_inc, sp2_inc)
        prob_share_site <- rep(0, (nsite + 1))
        if (prob == "hyper") {
            if (only_effects == FALSE) {
                all.probs <- phyper(0:min_inc, min_inc, nsite - 
                  min_inc, max_inc)

                prob_share_site[1] <- all.probs[1]
                prob_share_site[2:length(all.probs)] <- diff(all.probs)
            }
            else {

                ix1 <- (sp1_inc + sp2_inc) <= (nsite + 0:nsite) 
                ix2 <- 0:nsite <= min_inc
                prob_share_site[ix1 & ix2] <- 1
            }
        }
        if (prob == "comb") {
            if (only_effects == FALSE) {

                ix1 <- (sp1_inc + sp2_inc) <= (nsite + 0:nsite) 
                ix2 <- 0:nsite <= min_inc
                prob_share_site[ix1 & ix2] <- coprob(max_inc = max_inc, 
                        j = (0:nsite)[ix1 & ix2], min_inc = min_inc, nsite = nsite)

            }
            else {

                ix1 <- (sp1_inc + sp2_inc) <= (nsite + 0:nsite) 
                ix2 <- 0:nsite <= min_inc
                prob_share_site[ix1 & ix2] <- 1
            }
        }
        ix <- 0:nsite <= obs_cooccur[row, 3]
        p_lt <- sum(prob_share_site[ix])
        ix <- (0:nsite >= obs_cooccur[row, 3])
        p_gt <- sum(prob_share_site[ix])
        ix <- (0:nsite == obs_cooccur[row, 3])
        p_exactly_obs <- prob_share_site[max(which(ix))]
        p_lt <- round(p_lt, 5)
        p_gt <- round(p_gt, 5)
        p_exactly_obs <- round(p_exactly_obs, 5)
        prob_cooccur[row, 3] <- round(prob_cooccur[row, 3], 3)
        exp_cooccur[row, 3] <- round(exp_cooccur[row, 3], 1)
        #setTxtProgressBar(pb, nspp + row)
        return(c(sp1, sp2, sp1_inc, sp2_inc, obs_cooccur[row, 
            3], prob_cooccur[row, 3], exp_cooccur[row, 3], p_lt, 
            p_gt))       
        }))
    output <- as.data.frame(output)
    colnames(output) <- c("sp1", "sp2", "sp1_inc", "sp2_inc", 
        "obs_cooccur", "prob_cooccur", "exp_cooccur", "p_lt", 
        "p_gt")

    #close(pb)
    if (spp_names == TRUE) {
        sp1_name <- merge(x = data.frame(order = 1:length(output$sp1), 
            sp1 = output$sp1), y = spp_key, by.x = "sp1", by.y = "num", 
            all.x = T, sort = FALSE)
        sp2_name <- merge(x = data.frame(order = 1:length(output$sp2), 
            sp2 = output$sp2), y = spp_key, by.x = "sp2", by.y = "num", 
            all.x = T, sort = FALSE)
        output$sp1_name <- sp1_name[with(sp1_name, order(order)), 
            "spp"]
        output$sp2_name <- sp2_name[with(sp2_name, order(order)), 
            "spp"]
    }
    true_rand <- (nrow(output[(output$p_gt >= 0.05 & output$p_lt >= 
        0.05) & (abs(output$obs_cooccur - output$exp_cooccur) <= 
        (nsite * true_rand_classifier)), ]))
    output_list <- list(call = match.call(), results = output, 
        positive = nrow(output[output$p_gt < 0.05, ]), negative = nrow(output[output$p_lt < 
            0.05, ]), co_occurrences = (nrow(output[output$p_gt < 
            0.05 | output$p_lt < 0.05, ])), pairs = nrow(output), 
        random = true_rand, unclassifiable = nrow(output) - (true_rand + 
            nrow(output[output$p_gt < 0.05, ]) + nrow(output[output$p_lt < 
            0.05, ])), sites = nsite, species = nspp, percent_sig = (((nrow(output[output$p_gt < 
            0.05 | output$p_lt < 0.05, ])))/(nrow(output))) * 
            100, true_rand_classifier = true_rand_classifier)
    if (spp_names == TRUE) {
        output_list$spp_key <- spp_key
        output_list$spp.names = row.names(spp_site_mat)
    }
    else {
        output_list$spp.names = c(1:nrow(spp_site_mat))
    }
    if (thresh == TRUE) {
        output_list$omitted <- n_omitted
        output_list$pot_pairs <- n_pairs
    }
    class(output_list) <- "cooccur"
    if (only_effects == F) {
        output_list
    }
    else {
        effect.sizes(mod = output_list, standardized = eff_standard, 
            matrix = eff_matrix)
    }
}
