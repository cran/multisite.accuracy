multisite.accuracy <-
function (y, y.pred, site, estimate = c("auc", "bac", "cor", 
    "hr", "mse"), site.method = "covar", mixed = FALSE, min.n = 10, 
    ...) 
{
    if (!((is.vector(estimate) || is.factor(estimate)) && length(estimate) == 
        1 && estimate %in% c("auc", "bac", "cor", "hr", "mse"))) {
        stop("estimate must be \"auc\", \"bac\", \"cor\", \"hr\", or \"mse\"")
    }
    if (estimate %in% c("auc", "bac")) {
        if (!(is.vector(y) && all(y %in% 0:1))) {
            stop("for \"auc\" and \"bac\", y must be a binary vector")
        }
        n = length(y)
    }
    if (estimate %in% c("cor", "mse")) {
        if (!(is.vector(y) && is.numeric(y))) {
            stop("for \"cor\" and \"mse\", y must be a numeric vector")
        }
        n = length(y)
    }
    if (estimate == "hr") {
        if (!(class(y) == "Surv")) {
            stop("for \"hr\", y must be a Surv object")
        }
        n = nrow(y)
    }
    if (estimate %in% c("auc", "cor", "hr", "mse")) {
        if (!(is.vector(y.pred) && is.numeric(y.pred) && length(y.pred) == 
            n)) {
            stop("for \"auc\", \"cor\", hr\", and \"mse\", y.pred must be a numeric vector with the same length as y")
        }
    }
    if (estimate == "bac") {
        if (!(is.vector(y.pred) && all(y.pred %in% 0:1) && length(y.pred) == 
            n)) {
            stop("for \"bac\", y.pred must be a binary vector with the same length as y")
        }
    }
    if (!((is.vector(site.method) || is.factor(site.method)) && 
        length(site.method) == 1 && site.method %in% c("covar", 
        "meta", "none"))) {
        stop("site.method must be \"covar\", \"meta\" or \"none\"")
    }
    if (site.method != "none") {
        if (!((is.vector(site) || is.factor(site)) && length(site) == 
            n)) {
            stop("site must be a vector with the same length as y")
        }
        n.site = table(site)
        i_discard = which(site %in% names(n.site[which(n.site < 
            min.n)]))
        if (length(i_discard) > 0) {
            y = y[-i_discard]
            y.pred = y.pred[-i_discard]
            site = site[-i_discard]
        }
        if (site.method == "covar") {
            site = as.factor(site)
        }
    }
    if (!(is.numeric(min.n) && length(min.n) == 1 && min.n > 
        0)) {
        stop("min.n must be a positive number")
    }
    if (!(is.logical(mixed) && length(mixed) == 1)) {
        stop("mixed must be TRUE or FALSE")
    }
    if (estimate == "auc") {
        auc = switch(site.method, covar = {
            m = AROC.sp(formula.healthy = y.pred ~ site, group = "y", 
                tag.healthy = 0, data = data.frame(y, y.pred, 
                  site), ...)
            auc = unname(m$AUC[1])
            attr(auc, "method") = "AROC.sp"
            auc
        }, meta = {
            yi = vi = c()
            for (sitej in unique(site)) {
                i_sitej = which(site == sitej)
                mj = roc(y[i_sitej], y.pred[i_sitej], levels = 0:1, 
                  direction = "<")
                auc_cij = c(ci(mj))
                yi = c(yi, auc_cij[2])
                vi = c(vi, (max(diff(auc_cij))/qnorm(0.975))^2)
            }
            m = rma(yi = yi, vi = vi, ...)
            auc = c(m$beta)
            attr(auc, "method") = "roc-rma"
            auc
        }, none = {
            m = roc(y, y.pred, levels = 0:1, direction = "<")
            auc = c(m$auc)
            attr(auc, "method") = "roc"
            auc
        })
        return(data.frame(auc, site.method, auc.method = attr(auc, 
            "method")))
    }
    if (estimate == "bac") {
        for (group in 1:0) {
            y.pred.g = y.pred[which(y == group)]
            site.g = site[which(y == group)]
            se_sp = switch(site.method, covar = {
                contrasts(site.g) <- contr.sum
                m = logistf(y.pred.g == group ~ site.g, ...)
                beta0 = unname(m$coefficients[1])
                attr(beta0, "method") = "logistf"
                plogis(beta0)
            }, meta = {
                xi = ni = c()
                for (sitej in unique(site.g)) {
                  i_sitej = which(site.g == sitej)
                  xi = c(xi, sum(y.pred.g[i_sitej] == group))
                  ni = c(ni, length(i_sitej))
                }
                m = rma(xi = xi, ni = ni, measure = "PLO", ...)
                logit.se_sp = c(m$beta)
                attr(logit.se_sp, "method") = "rma/logit"
                plogis(logit.se_sp)
            }, none = mean(y.pred[which(y == group)] == group))
            if (group) {
                se = se_sp
            }
            else {
                sp = se_sp
            }
        }
        return(data.frame(bac = (se + sp)/2, se, sp, site.method, 
            se.method = ifelse(!is.null(attr(se, "method")), 
                attr(se, "method"), NA), sp.method = ifelse(!is.null(attr(sp, 
                "method")), attr(sp, "method"), NA), se.warning = ifelse(!is.null(attr(se, 
                "warning")), attr(se, "warning"), NA), sp.warning = ifelse(!is.null(attr(sp, 
                "warning")), attr(sp, "warning"), NA)))
    }
    if (estimate == "cor") {
        cor = switch(site.method, covar = {
            fe = function() {
                m = lm(y ~ y.pred + site, ...)
                t_df = c(summary(m)$coefficients[2, 3], m$df.residual)
                attr(t_df, "method") = "lm"
                t_df
            }
            me = function() {
                m = lmer(y ~ y.pred + (1 | site), control = lmerControl(check.conv.singular = .makeCC(action = "warning", 
                  tol = formals(isSingular)$tol)), ...)
                t_df = unname(summary(m)$coefficients[2, 4:3])
                attr(t_df, "method") = "lmer"
                t_df
            }
            me_error = function(msg) {
                t_df = fe()
                attr(t_df, "warning") = paste0("lmer returned an error/warning/message - ", 
                  msg$message)
                t_df
            }
            if (mixed) {
                t_df = tryCatch({
                  me()
                }, error = me_error, warning = me_error)
            } else {
                t_df = fe()
            }
            cor = 1/sqrt(1 + t_df[2]/t_df[1]^2)
            attr(cor, "method") = attr(t_df, "method")
            cor
        }, meta = {
            ri = ni = c()
            for (sitej in unique(site)) {
                i_sitej = which(site == sitej)
                ri = c(ri, cor(y.pred[i_sitej], y[i_sitej]))
                ni = c(ni, length(i_sitej))
            }
            m = rma(ri = ri, ni = ni, measure = "ZCOR", ...)
            z.cor = c(m$beta)
            attr(z.cor, "method") = "rma/fisher"
            tanh(z.cor)
        }, none = cor(y.pred, y))
        return(data.frame(cor, site.method, cor.method = ifelse(!is.null(attr(cor, 
            "method")), attr(cor, "method"), NA), cor.warning = ifelse(!is.null(attr(cor, 
            "warning")), attr(cor, "warning"), NA)))
    }
    if (estimate == "hr") {
        hr = switch(site.method, covar = {
            fe = function() {
                m = coxph(y ~ y.pred + site, ...)
                beta1 = unname(m$coefficients[1])
                attr(beta1, "method") = "coxph"
                beta1
            }
            me = function() {
                m = coxme(y ~ y.pred + (1 | site), ...)
                beta1 = unname(m$coefficients)
                attr(beta1, "method") = "coxme"
                beta1
            }
            me_error = function(msg) {
                beta1 = fe()
                attr(beta1, "warning") = paste0("coxme returned an error/warning/message - ", 
                  msg$message)
                beta1
            }
            if (mixed) {
                beta1 = tryCatch({
                  me()
                }, error = me_error, warning = me_error)
            } else {
                beta1 = fe()
            }
            exp(beta1)
        }, meta = {
            yi = vi = c()
            for (sitej in unique(site)) {
                i_sitej = which(site == sitej)
                mj = coxph(y[i_sitej] ~ y.pred[i_sitej])
                yi = c(yi, unname(mj$coefficients))
                vi = c(vi, c(mj$var))
            }
            m = rma(yi = yi, vi = vi, ...)
            log.hr = c(m$beta)
            attr(log.hr, "method") = "coxph-rma/log"
            exp(log.hr)
        }, none = {
            m = coxph(y ~ y.pred)
            beta1 = unname(m$coefficients)
            attr(beta1, "method") = "coxph"
            exp(beta1)
        })
        return(data.frame(hr, site.method, hr.method = attr(hr, 
            "method")))
    }
    if (estimate == "mse") {
        for (type in c("mean", "pred")) {
            mse.mean_pred = switch(site.method, covar = {
                fe = function() {
                  m = switch(type, mean = lm(y ~ site, ...), 
                    pred = lm(y ~ offset(y.pred) + site, ...))
                  e = unname(m$residuals)
                  attr(e, "method") = "lm"
                  e
                }
                me = function() {
                  m = switch(type, mean = lmer(y ~ (1 | site), 
                    control = lmerControl(check.conv.singular = .makeCC(action = "warning", 
                      tol = formals(isSingular)$tol)), ...), 
                    pred = lmer(y ~ offset(y.pred) + (1 | site), 
                      control = lmerControl(check.conv.singular = .makeCC(action = "warning", 
                        tol = formals(isSingular)$tol)), ...))
                  e = m@resp$y - m@resp$mu
                  attr(e, "method") = "lmer"
                  e
                }
                me_error = function(msg) {
                  e = fe()
                  attr(e, "warning") = paste0("lmer returned an error/warning/message - ", 
                    msg$message)
                  e
                }
                if (mixed) {
                  e = tryCatch({
                    me()
                  }, error = me_error, warning = me_error)
                } else {
                  e = fe()
                }
                mse.mean_pred = mean(e^2)
                attr(mse.mean_pred, "method") = attr(e, "method")
                attr(mse.mean_pred, "warning") = attr(e, "warning")
                mse.mean_pred
            }, meta = {
                log.se.mean_pred = c()
                yi = vi = c()
                for (sitej in unique(site)) {
                  i_sitej = which(site == sitej)
                  ej = y[i_sitej] - switch(type, mean = mean(y[i_sitej]), 
                    pred = y.pred[i_sitej])
                  msej = mean(ej^2)
                  vij = 1/(2 * (length(i_sitej) - 1))
                  yi = c(yi, log(sqrt(msej)) + vij)
                  vi = c(vi, vij)
                }
                m = rma(yi = yi, vi = vi, ...)
                log.se.mean_pred = c(m$beta)
                attr(log.se.mean_pred, "method") = "rma/nakagawa"
                exp(2 * log.se.mean_pred)
            }, none = switch(type, mean = mean((y - mean(y))^2), 
                pred = mean((y - y.pred)^2)))
            if (type == "mean") {
                mse.mean = mse.mean_pred
            }
            if (type == "pred") {
                mse.pred = mse.mean_pred
            }
        }
        return(data.frame(mse.mean, mse.pred, mse.pred_div_mse.mean = mse.pred/mse.mean, 
            site.method, mse.mean.method = ifelse(!is.null(attr(mse.mean, 
                "method")), attr(mse.mean, "method"), NA), mse.pred.method = ifelse(!is.null(attr(mse.pred, 
                "method")), attr(mse.pred, "method"), NA), mse.mean.warning = ifelse(!is.null(attr(mse.mean, 
                "warning")), attr(mse.mean, "warning"), NA), 
            mse.pred.warning = ifelse(!is.null(attr(mse.pred, 
                "warning")), attr(mse.pred, "warning"), NA)))
    }
}
