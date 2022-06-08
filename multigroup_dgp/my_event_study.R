# Estimate event-study coefficients using TWFE and proposed improvements.

library(staggered)

my_event_study = function(data, yname, idname, gname, tname, xformla = NULL, horizon = NULL, weights = NULL){
  
  # Setup ------------------------------------------------------------------------
  
  # Treat
  data$zz000treat = 1 * (data[[tname]] >= data[[gname]]) * (data[[gname]] > 0)
  data[is.na(data$zz000treat), "zz000treat"] = 0
  
  # Set g to zero if NA
  data[is.na(data[[gname]]), gname] = 0
  
  # Create event time
  data = data %>% dplyr::mutate(
    zz000event_time = dplyr::if_else(
      is.na(!!rlang::sym(gname)) | !!rlang::sym(gname) == 0,
      -Inf,
      as.numeric(!!rlang::sym(tname) - !!rlang::sym(gname))
    )
  )
  
  
  event_time = unique(data$zz000event_time)
  event_time = event_time[!is.na(event_time) & is.finite(event_time)]
  # All horizons
  if(is.null(horizon)) horizon = event_time
  
  # Format xformla for inclusion
  if(!is.null(xformla)) {
    xformla_null = paste0("0 + ", as.character(xformla)[[2]])
  } else {
    xformla_null = "0"
  }
  
  # TWFE -------------------------------------------------------------------------
  
  cli::cli_text("Estimating TWFE Model")
  
  tidy_twfe = NULL
  
  try({
    twfe_formula = stats::as.formula(glue::glue("{yname} ~ 1 + {xformla_null} + i(zz000event_time, ref = -1) | {idname} + {tname}"))
    est_twfe = fixest::feols(twfe_formula, data = data, warn = F, notes = F)
    
    tidy_twfe = broom::tidy(est_twfe) %>%
      dplyr::filter(stringr::str_detect(term, "zz000event_time::")) %>%
      dplyr::mutate(
        term = stringr::str_replace(term, "zz000event_time::", ""),
        term = as.numeric(term)
      ) %>%
      dplyr::select(term, estimate, std.error) %>%
      dplyr::mutate(estimator = "TWFE")
  })
  
  if(is.null(tidy_twfe)) cli::cli_warn("TWFE Failed")
  
  # did2s ------------------------------------------------------------------------
  
  cli::cli_text("Estimating using Gardner (2021)")
  
  tidy_did2s = NULL
  
  try({
    did2s_first_stage = stats::as.formula(glue::glue("~ 0 + {xformla_null} | {idname} + {tname}"))
    
    est_did2s = did2s::did2s(data, yname = yname, first_stage = did2s_first_stage, second_stage = ~i(zz000event_time, ref=-Inf), treatment = "zz000treat", cluster_var = idname, verbose = FALSE)
    
    tidy_did2s = broom::tidy(est_did2s) %>%
      dplyr::filter(stringr::str_detect(term, "zz000event_time::")) %>%
      dplyr::mutate(
        term = stringr::str_replace(term, "zz000event_time::", ""),
        term = as.numeric(term)
      ) %>%
      dplyr::select(term, estimate, std.error) %>%
      dplyr::mutate(estimator = "Gardner (2021)")
  })
  
  if(is.null(tidy_did2s)) cli::cli_warn("Gardner (2021) Failed")
  
  
  # did --------------------------------------------------------------------------
  
  cli::cli_text("Estimating using Callaway and Sant'Anna (2020)")
  
  tidy_did = NULL
  
  try({
    est_did =  did::aggte(
      did::att_gt(yname = yname, tname = tname, idname = idname, gname = gname, data = data)
        , type = "dynamic", na.rm = TRUE)
    
    tidy_did = data.frame("term" = est_did$egt, "estimate" = est_did$att.egt, "std.error"= est_did$se.egt, 
                     "estimator" = "Callaway and Sant'Anna (2020)") 
  })
  
  if(is.null(tidy_did)) cli::cli_warn("Callaway and Sant'Anna (2020) Failed")
  
  # sunab ------------------------------------------------------------------------

  # cli::cli_text("Estimating using Sun and Abraham (2020)")
  # 
  # tidy_sunab = NULL
  # 
  # try({
  #   # Format xformla for inclusion
  #   if(is.null(xformla)) {
  #     sunab_xformla = "1"
  #   } else {
  #     sunab_xformla = paste0("1 + ", as.character(xformla)[[2]])
  #   }
  # 
  #   sunab_formla = stats::as.formula(glue::glue("{yname} ~ {sunab_xformla} + sunab({gname}, {tname})"))
  # 
  #   est_sunab = fixest::feols(sunab_formla, data = data)
  # 
  #   tidy_sunab = broom::tidy(est_sunab) %>%
  #     filter(stringr::str_detect(term, glue::glue("{tname}::"))) %>%
  #     dplyr::mutate(
  #       term = stringr::str_replace(term, glue::glue("{tname}::"), ""),
  #       term = as.numeric(term)
  #     ) %>%
  #     dplyr::select(term, estimate, std.error) %>%
  #     dplyr::mutate(estimator = "Sun and Abraham (2020)")
  # })
  # 
  # if(is.null(tidy_sunab)) cli::cli_warn("Sun and Abraham (2020) Failed")

  # did_imputation ---------------------------------------------------------------

  cli::cli_text("Estimating using Borusyak, Jaravel, Spiess (2021)")

  tidy_impute = NULL

  try({
    impute_first_stage = stats::as.formula(glue::glue("~ {xformla_null} | {idname} + {tname}"))

    tidy_impute = did2s::did_imputation(data,
                                        yname = yname, gname = gname, tname = tname, idname = idname,
                                        first_stage = impute_first_stage, horizon = TRUE, pretrends = TRUE) %>%
      dplyr::select(term, estimate, std.error) %>%
      dplyr::mutate(estimator = "Borusyak, Jaravel, Spiess (2021)", term = as.numeric(term))
  })

  if(is.null(tidy_impute)) cli::cli_warn("Borusyak, Jaravel, Spiess (2021) Failed")

  # staggered --------------------------------------------------------------------

  # cli::cli_text("Estimatng using Roth and Sant'Anna (2021)")
  # 
  # tidy_staggered = NULL
  # 
  # try({
  #   tidy_staggered = staggered::staggered(
  #     data %>%
  #       dplyr::mutate(!!rlang::sym(gname) := dplyr::if_else(!!rlang::sym(gname) == 0, Inf, !!rlang::sym(gname))),
  #     i = idname, t = tname, g = gname, y = yname, estimand = "eventstudy",
  #     eventTime = event_time[is.finite(event_time) & event_time != -1]
  #   ) %>%
  #     dplyr::select(term = eventTime, estimate, std.error = se) %>%
  #     dplyr::mutate(estimator = "Roth and Sant'Anna (2021)")
  # })
  # 
  # if(is.null(tidy_staggered)) cli::cli_warn("Roth and Sant'Anna (2021) Failed")



  
  # Bind results together --------------------------------------------------------
  
  out = dplyr::bind_rows(tidy_twfe, tidy_did2s, tidy_did, #tidy_sunab, 
                         #tidy_staggered,
                         tidy_impute
                         )
  
  return(out)
  
}