# Input : give Class and Prefix/list of prefix
# Return list of corresponding columns
get_col_names <- function(class_dict, drugsets_dict, class = "c") {
  if (class == "c") {
   col_names = unique(class_dict$drugclass_c)
  } else if (class == "b") {
    col_names = unique(class_dict$drugclass_b)
    col_names <- substr(col_names,1,20)
  } else if (class == "a") {
    col_names = unique(class_dict$drugclass_a)
    col_names <- substr(col_names,1,20)
  } else if (class == "sets") {
    col_names = unique(drugsets_dict$drugset_alias)
    col_names <- substr(col_names,1,20)
  }  else {
    stop("Invalid class name. Should be either 'a', 'b', 'c' or 'sets'")
  }
  col_names <- col_names[col_names != ""]
  return(col_names)
}

# Return the column names corresponding to a drugclass/drugset and a list of prefix
get_chemo_column <- function(class_dict, drugsets_dict, class = "c", prefix = "pct", verbose = F) {
  # Get drugclasses names
  col_names <- get_col_names(class_dict = class_dict, drugsets_dict = drugsets_dict, class = class)  
  output <- c()
  # Get prefix
  if (length(prefix) >= 1) {
    for (pref in prefix) {
      if (verbose) {print(pref)}
      output <- output %>% c(paste0(pref, "_", col_names))
    } 
  }
  else {
    stop("Invalid prefix argument. Should either be a string or a list of string")
  }
  
  return(output)  
}

#get_chemo_column("c", c("adm_pmonth", "dose_pdamin", "length"))

# Return a list of string corresponding to interactions terms in a formula
get_chemo_column_interactions <- function(class = "c", interaction = NULL, interaction_type = ":") {
  # 'Interat'interaction' should be a list of arrays, each containing the desired group of features
  col_names <- get_col_names(class = class)  
  output <- c()
  if (length(interaction) > 0) {
    if (typeof(interaction) != "list") {stop("interaction must be a list")}
    for (inter_list in interaction) {
      stopifnot(length(inter_list) >= 2)
      for (col in col_names) {
      inter_term <- ""
        for (inter_pref in inter_list) {
          inter_term <- paste0(inter_term, interaction_type, inter_pref, "_", col)
        }
      # Supress the last ":"
      inter_term <- substr(inter_term, 2, nchar(inter_term))
      output <- output %>% c(inter_term)
      }
    }
    return(output)
  } else {
    return(NULL)
  }
}
  
#get_chemo_column_interactions(interaction = list(c("lenght", "pct"), c("ds_padm", "ind")))

remove_prefix <- function(col_names, prefix) {
  #' Remove for the list "col_names" every value starting with "prefix"
  for (pref in prefix) {
    col_names <- col_names[!(col_names %>%grepl(pattern = pref))]
  }
  return(col_names)
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
        filter(n >= threshold) %$%
        term
    
    excluded = setdiff(terms, terms_filtered)
    
    if (length(excluded) > 0 & verbose) {
        print(paste0(excluded, ' were excluded due to insufficient data', collapse = ', '))
    }
    
    return(terms_filtered)
    
}

# helper functions
term_count = function(x) {
    if (is_indicator(x)) {
        sum(x)
    } else {
        sum(!is.na(x))
    }
}

# check if vector of data is a binary indicator
is_indicator = function(x) {
    
    values = x %>% as.character() %>% unique
    
    return(
        setequal(values, c('0', '1')) | setequal(values, c('0')) | setequal(values, c('1'))
    )
}



