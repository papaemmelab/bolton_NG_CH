# Toolbox.R

require(reshape2)

signif.num <- function(x, ns = FALSE) {
    if (ns) {
        symbols = c("***", "**", "*", "ns")
    } else {
        symbols = c("***", "**", "*", "")
    }
    
    symnum(unlist(x), corr = FALSE, na = FALSE, legend = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
           symbols = symbols)
}

median0 = function(x) {
  if (length(x) == 0) {
    return(NA)
  } else {
    return(median(x))
  }
}

count_level = function(D, x) {
    if (is_continuous(D[[x]])) {
        R = data.frame(n_level = nrow(D)) %>% mutate(variable = x, levels = x)
    } else {
        R = D %>% 
          filter(!is.na(get(x))) %>%
          count(get(x)) %>%
          mutate(variable = x) %>%
          dplyr::rename(levels = 'get(x)', n_level = n)
    }
    R %>% select(variable, levels, n_level)
}

get_xs = function(model) {
    vars = model$terms %>% attr("variables") %>% as.character %>% .[-1]
    xs = vars[-1]
    return(xs)
}

is_integer = function(x) {
    if (all(is.na(as.numeric(x)))) {
        return(FALSE)
    } else {
        unique(as.numeric(x)) %>% .[!is.na(.)] %>% {. %% 1 == 0} %>% all
    }
}

is_continuous = function(x) {
    is.numeric(x) & (!is_integer(x))
}

# Filter drugs that do not have sufficient counts 
filter_terms = function(M_wide, terms, threshold = 10, verbose = T) {

    terms_filtered = M_wide %>% 
        summarise_at(
            terms,
            term_count
        ) %>% t %>%
        as.data.frame %>%
        tibble::rownames_to_column('term') %>%
        dplyr::rename(n = V1) %>%
        dplyr::filter(n >= threshold) %$%
        term
    
    excluded = setdiff(terms, terms_filtered)
    
    if (length(excluded) > 0 & verbose) {
        print(paste0(c(excluded, ' were excluded due to insufficient data'), collapse = ', '))
    }
    
    return(terms_filtered)
    
}

z_score = function(x) {
  return((x - mean(x)) / sd(x))
}

# helper functions
term_count = function(x) {
    if (is_indicator(x)) {
        sum(as.integer(as.character(x[!is.na(x)])))
    } else if (suppressWarnings(is_integer(x))) {
        sum(as.integer(x)[!is.na(x)] != 0)
    } else {
        sum(!is.na(x))
    }
}

# check if vector of data is a binary indicator
is_indicator = function(x) {

    x = x[!is.na(x)]
    
    values = x %>% as.character() %>% unique
    
    return(
        setequal(values, c('0', '1')) | setequal(values, c('0')) | setequal(values, c('1'))
    )
}

multinom_summary = function(model) {

#     confint(model) %>% as.data.frame
    
    p_vals = t(summary(model)$coefficients / summary(model)$standard.errors) %>%
        abs() %>% pnorm(lower.tail = F) %>% (function(x) {x*2}) %>%
        melt(
            varnames = c('coefficient', 'category'),
            value.name = 'p_val'
        ) %>%
        filter(coefficient != '(Intercept)')

    ORs = model %>% coef %>% t %>% exp %>%
        melt(
            varnames = c('coefficient', 'category'),
            value.name = 'OR'
        ) %>%
        filter(coefficient != '(Intercept)')
    
    summary = full_join(ORs, p_vals, by = c("coefficient", "category"))
    
    summary$signif = pval_to_signif_code(summary$p_val)
    
    summary = summary %>% select(category, coefficient, OR, p_val, signif)
    
    return(summary)
}

pval_to_signif_code <- function(p_vals) {
  return(ifelse(p_vals < 0.001, "***",
                ifelse(p_vals < 0.01, "**",
                       ifelse(p_vals < 0.05, "*",
                              ifelse(p_vals < 0.1, ".", "")))))
}

# Summary table of a regression model
# summar = function(logit) {
#     
#   summary = as.data.frame.matrix(
#     cbind(
#         Odds_Ratio=exp(coef(logit)),
#         SE = coef(summary(logit))[,2], exp(confint.default(logit)),
#         p_vals = coef(summary(logit))[,4])
#   )
#     
#   summary$Signif = pval_to_signif_code(summary$p_vals)
#   summary = summary[row.names(summary) != "(Intercept)", ]
#     
#     
#     
#   summary = summary %>% tibble::rownames_to_column(var = 'Predicator') %>% 
#     mutate_at(.vars = c('Odds_Ratio', 'SE', '2.5 %', '97.5 %', 'p_vals'), .funs = as.numeric) %>%
#     mutate_if(.predicate = is.numeric, .funs = function(x){signif(x, 2)})
# 
#   return(summary)
# }

ord_summary <- function(ord_logit, confint = TRUE) {
  # Deleting incept rows
  intercept_names <- names(summary(ord_logit)$zeta)
  coefs <- as.data.frame.matrix(coef(summary(ord_logit)))
  coefs <- coefs[!row.names(coefs) %in% intercept_names,]
  
  if (confint) {
    coef_table <- as.data.frame.matrix(cbind(Odds_Ratio=exp(coefs[,1]), 
                                             SE = coefs[,2], 
                                             exp(confint.default(ord_logit)), 
                                             p_vals = coefs[,3] %>%abs() %>% pnorm(lower.tail = FALSE) * 2))
  } else {
    coef_table <- as.data.frame.matrix(cbind(Odds_Ratio=exp(coefs[,1, drop = FALSE]), 
                                             SE = coefs[,2], 
                                             p_vals = coefs[,3] %>%abs() %>% pnorm(lower.tail = FALSE) * 2))
  }
  coef_table$Signif <- pval_to_signif_code(coef_table$p_vals)
  return(coef_table)
}

get_formula <- function(logit, clean = T) {
  # GLM formula are accessed by "formula", while Polr's are accessed by "terms"
  formula = ifelse(!is.null(logit$formula), Reduce(paste, deparse(logit$formula)), Reduce(paste, deparse(logit$terms)))
  if (clean) {
    return(clean_formula(formula))
  } else {
    return(formula)
  }
}

clean_formula <- function(formula) {
  formula <- formula %>% gsub(pattern = "\\)", replacement = "") %>% 
    gsub(pattern = "get\\(", replacement = "") %>%
    gsub(pattern = "factor\\(", replacement = "")
}


# Plotting tools
plot_summary <- function (summary, logit = NULL, title = "Regression Summary", digits = 2, add_ref = TRUE) {
  #' Constructs an array showing the summary of a regression model
  #' 'summary' should be an instance of the function summary(logit), where logit is a regression model
  #' Entering a logit argument displays the model's formula below the caption
  if (!is.null(logit)) {
    # GLM formula are accessed by "formula", while Polr's are accessed by "terms"
    formula = ifelse(!is.null(logit$formula), Reduce(paste, deparse(logit$formula)), Reduce(paste, deparse(logit$terms)))
    title = paste0(title, "<br>", formula)
  }
  
  ref_rows <- c()
  if (add_ref) {
    for (f in factor_var) {
      row_labels <- summary %>% row.names() %>% gsub(pattern = "factor\\(", replacement = "") %>% 
        gsub(pattern = "\\)", replacement = "") %>% gsub(pattern = ".*:.*", replacement = "")
      
      ##/A/ Clean all row labels
      row.names(summary) <- row_labels
      
      index <- which(row_labels %>% grepl(pattern = paste0("^", f, "[^_]")) == TRUE)
      if (length(index) > 0) {
        ref <- paste0(f, factor_dict[[f]])
        m_index <- min(index)
        summary <- rbind(summary[1:m_index-1, ], c(1, 0, 1, 1, 0, 0), summary[m_index:nrow(summary), ])
        # Apply verbose name to reference group
        ref <- ifelse(ref %in% names(factor_label_dict), factor_label_dict[[ref]], ref)
        row.names(summary)[m_index] <- paste0(ref, "  -  Ref")
        ref_rows <- cbind(ref_rows, m_index)
      }
    }
  }
  #digit_vector <- rep(digits, nrow(summary))
  hidden_row_spec <- ref_rows
  row.names(summary) <- lapply(row.names(summary), FUN = function(x) ifelse(x %in% names(factor_label_dict), 
                                                                            factor_label_dict[[x]], x))
  out_table <- kable(summary, caption = title, digits = digits) %>%
  #  kable_styling("striped", full_width = F) %>%
    kable_styling(full_width = F) %>%
    kable_styling(font_size = 16) %>%
    column_spec(c(2), bold = T) %>% 
    row_spec(c(1:nrow(summary)), color = "black") %>%
    column_spec(c(1:ncol(summary) + 1), color = "black")
  
  if (length(hidden_row_spec) > 0) {
    return(
      out_table %>% row_spec(hidden_row_spec, color = "black") %>%
        column_spec(c(1:2), color = "black")
    )
  }
  else{
    return(out_table %>% column_spec(c(1), color = "black"))
  }
}

plot_forest <- function(class_results, x = "class", y = "estimate", ymin = "conf.low", ymax = "conf.high",
                       label = "p.label", limits = NULL, breaks = waiver(), title = "", col = NULL, fill = NULL,
                        dodge_width = 0.8, outer_limit_arrows = FALSE, ps=3, eb_w=0.4, eb_s=0.4, or_s=4, OR=T, yinter = 1,
                        nudge = 0, base_plot = geom_blank(), bar_col = 'black') { 
    # suppressWarnings(ggthemr::ggthemr("fresh"))
  #' Forest plot of the coefficients in 'class_results'
  output_plot <- ggplot(class_results, aes_string(x = x, y = y, ymin = ymin, ymax = ymax, col = col, fill = fill)) + 
    base_plot +
    geom_hline(yintercept = yinter, color = "gray", linetype = "solid") +
    geom_errorbar(position = position_dodge(width = dodge_width), width = eb_w, size=eb_s) + 
    geom_point(position = position_dodge(width = dodge_width), size=ps) +
    geom_text(aes_string(label = label, vjust = nudge), position = position_dodge(width = dodge_width), color = 'black', size = or_s, alpha = .9) +
    coord_flip() +
    ggtitle(title)
  
  if (OR==T) {
  output_plot <-  output_plot +   
    scale_y_log10(limits = limits, breaks = breaks)
  }
  
  if (outer_limit_arrows) {
    stopifnot(length(limits) > 0)
    # Check for errorbar values outside 
    class_results[, "ymin_out"] <- ifelse(class_results[, ymin] < limits[1], limits[1], NA)
    class_results[, "linestyle_min"] <- ifelse(class_results[, ymin] < limits[1], "b", NA)
    class_results[, "ymax_out"] <- ifelse(class_results[, ymax] > limits[2], limits[2], NA)
    class_results[, "linestyle_max"] <- ifelse(class_results[, ymax] > limits[2], "b", NA)
    
    output_plot <- output_plot + geom_linerange(data = class_results, aes_string(x = x, ymin = "ymin_out", ymax = ymax), linetype = 3, position = position_dodge(width = dodge_width))
    output_plot <- output_plot + geom_linerange(data = class_results, aes_string(x = x, ymin = ymin, ymax = "ymax_out"), linetype = 3, position = position_dodge(width = dodge_width))
    output_plot <- output_plot + geom_linerange(data = class_results, aes_string(x = x, ymin = "ymin_out", ymax = "ymax_out"), linetype = 3, position = position_dodge(width = dodge_width))
    
    #output_plot <- output_plot + geom_segment(data = class_results, aes_string(x = x, xend = x, y = ymax, yend = paste0("ymin_out")), 
    #                                          arrow = arrow(length = unit(.2, "cm")), position = position_dodge(width = dodge_width), lineend = "round")
    #output_plot <- output_plot + geom_segment(data = class_results, aes_string(x = x, xend = x, y = ymin, yend = paste0("ymax_out")), 
    #                                          arrow = arrow(length = unit(.2, "cm")), position = position_dodge(width = dodge_width), lineend = "round")
  }
  return(output_plot)
}

plot_table_forest <- function(class_results, omit_columns = c(), title= "", digits = 2, 
                              index_col = c(1), signif_only = FALSE, omit_columns_base = c("group", "xpos", "xmin", "xmax", "p.label")) {
  #' Returns a summary of the coefficients in 'class_results'
  if (signif_only) {
    class_results <- class_results %>% filter(p.value < 0.05)
  }
  omit_columns <- c(omit_columns, omit_columns_base)
  class_results <- class_results[,-which(names(class_results) %in% omit_columns)]
  reorder_columns <- c(index_col, (1:length(class_results))[!1:length(class_results) %in% index_col])
  class_results <- class_results[, reorder_columns]
  output <- kable(class_results, caption = title, digits = digits) %>%
    kable_styling("striped", full_width = F)  %>% 
    column_spec(c(1:length(index_col)), bold = T)
  
  return(output)
}

summar = function(model, truncate_p = TRUE, p_het = NULL, CI = TRUE, signif = FALSE, n = FALSE, latex = TRUE, ref = TRUE, sep = F, format = TRUE) {
    
    if ('geese' %in% class(model)) {
      summary = ordgee_summary(model) %>%
        select(OR = Odds.Ratio, conf.low = OR.lower, conf.high = OR.upper, pval = p)
    } else if ('polr' %in% class(model)) {
      summary = ord_summary(model) 
    } else if ('geeglm' %in% class(model)) {
      summary = as.data.frame.matrix(
          cbind(
              OR=exp(coef(model)),
              exp(broom::confint_tidy(model)),
              pval = coef(summary(model))[,4])
        ) %>%
        tibble::rownames_to_column('levels')
    } else {
        summary = sjPlot::get_model_data(model=model,type="est") %>%
          select(levels = term, OR = estimate, conf.low, conf.high, pval = p.value, p.stars) %>%
          mutate(levels = as.character(levels))
    }
        
    # rounding and filter terms
    summary = summary %>%
        filter(!str_detect(levels, 'Inter:|Intercept')) %>%
        mutate_at(.vars = c('OR', 'conf.low', 'conf.high', 'pval'), .funs = as.numeric) %>%
        mutate_if(.predicate = is.numeric, .funs = function(x){signif(x, 2)})
    
    # extract the formula
    if (!is.null(model$formula)) {
        formula = model$formula
    } else {
        formula = model$call[[2]]
    }
    
    formula = formula %>% deparse() %>% as.character() %>% Reduce(paste, x = .)
    
    # add factor column
    factors = formula %>%
        str_remove_all(' |"|\\*') %>%
        strsplit('~') %>% 
        unlist() %>%
        .[2] %>%
        strsplit('\\+') %>% 
        unlist()
        
    if (!is.null(p_het)) {
        
        factors = factors %>% str_remove(p_het)
        
        summary = summary %>%
            filter(str_detect(levels, 'therapy_binarytreated:|:therapy_binarytreated')) %>%
            mutate(levels = str_remove(levels, 'therapy_binarytreated:|:therapy_binarytreated')) 
    }
    
    summary = summary %>% rowwise() %>% 
      mutate(
        levels = ifelse(is_indicator(model$data[[levels]]) & !str_detect(levels, '1'), paste0(levels, '1'), levels)
      ) %>%
      mutate(
          variable = levels %>% str_detect(factors) %>% which %>% tail(n = 1) %>% factors[.],
          levels = ifelse(
              levels %in% factors,
              levels,
              levels %>% str_remove(pattern = paste0(factors, collapse = '|'))
          )
      ) %>%
      ungroup() 
  
    # concat OR and CI interval
    summary = summary %>%
        mutate(CI = paste0(conf.low, '-', conf.high))

    if (!CI) {
        summary = summary %>% mutate(OR = paste(OR, CI))
    }
    
    if (n) {
        # add n column
        summary = lapply(
                get_xs(model),
                function(x) { count_level(model$data, x) }
            ) %>%
            Reduce(f = rbind, x = .) %>%
            right_join(summary, by = c('levels', 'variable'))
        summary
    }

    # get reference levels
    if (ref) {
      refs = lapply(model$xlevels, function(x){x[1]}) %>% {.[. != '0']}

      summary = summary %>% rowwise() %>%
        mutate(ref = ifelse(variable %in% names(refs), unlist(unname(refs[variable])), '')) %>%
        ungroup()
    }
    
    if (truncate_p) {
        summary = summary %>% 
            mutate(
                pval = ifelse(pval < 0.000001, '<1e-06', pval)
            )
    }

    summary = summary %>% mutate(levels = ifelse(as.character(levels) == '1', variable, levels))

    # clean up formatting
    if (format) {
        summary = summary %>%
            # mutate(levels = ifelse(as.character(levels) == '1', variable, levels)) %>%
            mutate(
                variable = format_variable(variable),
                levels = format_variable(levels)
            )
        if (ref) {
            summary = summary %>% mutate(ref = format_variable(ref))
        }
    }

    # paste reference levels
    if (ref) {
        summary = summary %>% rowwise() %>%
            mutate(variable = ifelse(ref != '', paste0(variable, ' (', ref, ')'), variable)) %>%
            ungroup()
    }

    if (n) {
      summary = summary %>% mutate(levels = paste0(levels, ' (', n_level, ')'))
    }

    # arrange columns
    if (CI) {
        final_cols = c('variable', 'levels', 'OR', 'CI', 'pval')
    } else {
        final_cols = c('variable', 'levels', 'OR', 'pval')
    }
    
    if (signif) {
        final_cols = c(final_cols, 'p.stars')
    }

    summary = summary %>% select(final_cols)
    
    if (latex) {
      summary = summary %>%
        dplyr::rename(`95\\% CI` = `CI`) %>%
        mutate(pval = kableExtra::cell_spec(pval, "latex", bold = ifelse(as.numeric(pval) < 0.05 | pval == '<1e-06', T, F)))
    }

    if (sep) {
      summary = summary %>% tidyr::separate(CI, c('lower', 'upper'), '-') %>%
        mutate(OR = as.numeric(OR), upper = as.numeric(upper), lower = as.numeric(lower))
    }

    return(summary)
}

ordgee_summary <- function(fit, alpha=.05, dig=3, p.dig=4){
    zq <- qnorm(1-alpha/2)
    estimate <- fit$beta
    lower    <- fit$beta - zq*summary(fit)$mean[,"san.se"]
    upper    <- fit$beta + zq*summary(fit)$mean[,"san.se"]
    
    robust.se <- round(summary(fit)$mean[,"san.se"], dig)
    p <- summary(fit)$mean[,"p"]
    p <- signif(p, digits=p.dig)
    if(all(fit$model[1:2]==c("logit","binomial"))){
        Odds.Ratio <- round(exp(estimate), dig)
        OR.lower   <- round(exp(lower), dig)
        OR.upper   <- round(exp(upper), dig)
        estimate   <- round(estimate, dig)
        return(as.data.frame(cbind(estimate,robust.se,Odds.Ratio,OR.lower,OR.upper,p)))
    } else if(all(fit$model[1:2]==c("log","poisson"))){
        #assumes all counts are taken over same lengths of time 
        #or that an offset term has been specified, as usual.
        IRR       <- round(exp(estimate), dig) #incidence rate ratio
        IRR.lower <- round(exp(lower), dig)
        IRR.upper <- round(exp(upper), dig)
        estimate  <- round(estimate, dig)
        return(as.data.frame(cbind(estimate,robust.se,IRR,IRR.lower,IRR.upper,p)))
    } else{
        estimate <- round(estimate, dig)
        lower    <- round(lower, dig)
        upper    <- round(upper, dig)
        return(as.data.frame(cbind(estimate,robust.se,lower,upper,p)))
    }
}

ord_summary = function(model) {
    intercept_names = names(summary(model)$zeta)
    coefs = as.data.frame.matrix(coef(summary(model)))
    coefs = coefs[!row.names(coefs) %in% intercept_names,]

    as.data.frame(
        cbind(OR = exp(coefs[,1]), 
        as.data.frame(exp(confint.default(model))) %>% dplyr::rename(conf.low = `2.5 %`, conf.high = `97.5 %`), 
        pval = coefs[,3] %>%abs() %>% pnorm(lower.tail = FALSE) * 2)
    ) %>% tibble::rownames_to_column('levels')
}

format_variable = function(vars) {
  as.character(vars) %>%
    str_remove_all(regex('binary|scaled|bin|verbose|_b$| b$|^g_|^ind_|^pct_|_c$', ignore_case = T)) %>%
    str_replace_all('_', " ") %>%
    str_replace_all('-$', " ") %>%
    str_replace_all('^M$', "Male") %>%
    str_replace_all('^F$', "Female") %>%
    str_replace_all('SystemTumorType', 'Tumor Type') %>%
    str_trim() %>%
    tolower() %>%
    tools::toTitleCase() %>%
    cap_susbet(cap_set)
}

cap_set = c('DNA', '^HR$', 'HRD', 'MMR', 'CHEK2', 'CH', 'PD', 'II',
 'DNMT3A', 'TET2', 'PPM1D', 'ASXL1', 'TP53', 'ATM', 'SF3B1', 'JAK2', 'SRSF2', '^VAF',
 '^ANC$', '^ALC$', '^AMC$', '^HGB$', '^MCV$', '^PLT$', '^WBC$', '^RDW$', 'XRT', 'IDH')

# capitalize a subset of words in a vector of string
cap_susbet = function(x, words) {
    for (word in words) {
      x = str_replace(x, regex(word, ignore_case = T), str_remove_all(word, regex("\\^|\\$| ")))
    }
    return(x)
}

# display_kable = function(kable, file = 'test.png')  {
    
#     kable %>% save_kable(file = file)
    
#     if (tools::file_ext(file) == 'pdf') {
#         display_svg(file = pdf_to_svg(file))
#     } else {
#         display_png(file = file)
#     }
# }

pdf_to_svg = function(pdf_file) {
    
    svg_file = paste0(tools::file_path_sans_ext(test_table), '.svg')
    
    cmd = paste('pdf2svg', pdf_file, svg_file)
    
    msg = system(cmd, intern = TRUE)
        
    return(svg_file)
}


table_regression = function(summary, width = 80, caption = '')  {
        
    table = summary %>%
        kable(format = "latex", booktabs = T) %>%
        kable_styling(
            latex_options = c("striped", "hold_position"),
            full_width = T,
            font_size = 10
        ) %>%
        column_spec(2:50, width = width) %>%
        column_spec(1, width = 20) %>%
        collapse_rows(columns = 1:2, row_group_label_position = 'stack') 
    
    return(table)
}

crosstab = function(D, x, y, fraction = FALSE) {
        
    ylevels = unique(D[[y]]) %>% .[!is.na(.)] %>% as.character
    
    if (tolower(x) == 'total') {
        D %>% dcast(
            paste0('. ~ ', y), 
            fun.aggregate = length,
            value.var = y
        ) %>%
        dplyr::rename('levels' = '.') %>%
        mutate(levels = 'total') %>%
        mutate(variable = 'total') %>%
        select(variable, levels, ylevels)
    } else if (tolower(y) == 'total') {
        D %>% count(get(x)) %>% 
        dplyr::rename('levels' = 'get(x)', 'total' = 'n') %>%
        mutate(variable = gsub('_', '-', x)) %>%
        select(variable, levels, total)
    } else if (is_continuous(D[[y]]) & (is_continuous(D[[x]]) | class(D[[x]]) == 'integer')) {
        NA
    } else if ((is_continuous(D[[x]]) | class(D[[x]]) == 'integer') & !is_indicator(D[[x]])) {
        # x is continuous
        D %>% 
        filter(!is.na(get(x))) %>%
        dcast(
            paste0('. ~ ', y), 
            fun.aggregate = function(x){
                paste0(signif(mean(x), 2), ' (', signif(sd(x), 2), ')')
            },
            value.var = x
        ) %>%
        dplyr::rename('levels' = '.') %>%
        mutate(levels = 'mean (SD)') %>%
        mutate(variable = gsub('_', '-', x)) %>%
        select(variable, levels, ylevels)
    } else if (is_continuous(D[[y]])) {
        # y is continous
        D %>% dcast(
            paste0('. ~ ', x), 
            fun.aggregate = function(x){
                paste0(signif(mean(x), 2), ' (', signif(sd(x), 2), ')')
            },
            value.var = y
        ) %>% .[,-1] %>% t %>% 
        as.data.frame %>%
        tibble::rownames_to_column(var = 'levels') %>%
        setNames(c('levels', 'mean (SD)')) %>%
        mutate(variable = gsub('_', '-', x)) %>%
        select('variable', 'levels', 'mean (SD)')
    } else {
        res = table(D[[x]], D[[y]]) %>%
            as.data.frame.matrix %>% as.data.frame %>%
            tibble::rownames_to_column(var = 'levels') %>%
            mutate(variable = gsub('_', '-', x)) %>%
            select(variable, levels, ylevels)

        if (fraction) {
            res = res %>%
                group_by(variable) %>%
                mutate_at(ylevels, function(x){paste0(x, ' (', signif(x * 100/sum(x), 2), '%', ')')}) %>%
                ungroup()
        }

        return(res)
    }
}

tab_response = function(D, xs, y, fraction = FALSE) {    
    lapply(
        xs,
        function(x) {crosstab(D, x, y, fraction)}
    ) %>% 
    Reduce(rbind, .)
}