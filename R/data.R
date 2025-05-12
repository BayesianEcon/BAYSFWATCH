
#' Inflation and Unemployment Dataset
#'
#' This dataset contains monthly macroeconomic variables for 11 EU countries
#' (Austria, Belgium, Estonia, Finland, France, Germany, Greece, Ireland, Italy,
#' The Netherlands, Portugal, and Spain) from January 2002 to October 2022. The
#' variables include the Industrial Production Index, Price Index, and
#' Unemployment Rate, sourced from FRED and EUROSTAT. Series are adjusted for
#' seasonality and working days, and stationarity is enforced via differencing.
#'
#' Outliers are detected in key periods such as the COVID-19 pandemic (Mar–Aug 2020, Nov 2020),
#' the Ukrainian conflict (Mar–Apr 2022), and the inflation rise (Sep–Oct 2022).
#'
#' @format An array of dimension 3 x 11 x 250
#' @source Canova and Ciccarelli (2009) DOI: 10.1111/j.1468-2354.2009.00554.x
"EUdata"

#' Trade Network Dataset
#'
#' This dataset comprises 22 annual trade networks among 27 countries from 1995
#' to 2017, derived from IMF reports and other official sources (COMTRADE,
#' EUROSTAT). Trade values are transformed to constant USD (base year 2010) using
#' the US GDP PPP deflator. The selected countries are large economies, filtered
#' for data completeness.
#'
#' Notable outliers are detected in 2015 and 2017, while 2016 is not flagged as
#' anomalous.
#'
#' @format An array of dimension 27 x 27 x 22
#' @source Rose (2004) DOI: 10.1257/000282804322970724
"Tradedata"

#' Volatility Network Dataset
#'
#' This dataset consists of 145 weekly volatility networks (Monday–Friday) from
#' January 2016 to September 2020, constructed using pairwise Granger causality
#' for 50 firms from Germany, France, and Italy. Firms span 11 GICS sectors,
#' including Financials, IT, Industrials, Health Care, and others.
#'
#' Outliers capture volatility regime changes before and during the COVID-19
#' crisis, particularly transitions such as from "low" to "high" and "moderate"
#' volatility levels.
#'
#' @format An array of dimension 50 x 50 x 145
#' @source Billio, Casarin, and Iacopini (2022) DOI: 10.1080/07350015.2022.2032721
"Volatilitydata"
